import os
import subprocess
import tempfile
import sys
from numpy import quantile
from collections import defaultdict, Counter
from tqdm import tqdm
import logging
from .peeking_sort import interleave_streams, buffered_aggregator, filter_intersection
logger = logging.getLogger('Windows')
logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO)

class RegionOverlapComparitor:
    '''
    Why does this class exist?

    This class is used to wrap a region and change the comparison behavior of the region.
    When wrapped, a region is equal to another region if the two regions overlap.
    '''
    @classmethod
    def unwrap(cls, wrapped):
        return wrapped.unwrap()

    def __init__(self, chrom, start, end, *args):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.args = args

    def __gt__(self, other):
        return self.chrom > other.chrom or \
            (self.chrom == other.chrom and self.start > other.end)
    
    def __eq__(self, other):
        '''
        self  ____
        other   ____

        self      ____
        other  ____
        '''
        return self.chrom == other.chrom and \
                    (self.end >= other.start and self.start <= other.end)
    
    def __ge__(self, other):
        return (self > other) or (self == other)
    
    def unwrap(self):
        return self.chrom, self.start, self.end, *self.args


def check_regions_file(regions_file):

    encountered_idx = defaultdict(lambda : False)

    with open(regions_file, 'r') as f:

        for i, line in enumerate(f):
            
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            assert len(cols) >= 12, \
                f'Expected 12 or more columns (in BED12 format) in {regions_file}, with the fourth column being an integer ID.\n'
            
            try:
                _, start, end, idx = cols[:4]
                idx = int(idx); start = int(start); end = int(end)
            except TypeError:
                raise TypeError(
                    f'Count not parse line {i} in {regions_file}: {line}.\n'
                    'Make sure the regions file is tab-delimited with columns chr<str>, start<int>, end<int>, idx<int>.\n'
                )
            
            encountered_idx[idx] = True

    assert len(encountered_idx) > 0, \
        f'Expected regions file to have at least one region.'

    largest_bin = max(encountered_idx.keys())
    assert all([encountered_idx[i] for i in range(largest_bin + 1)]), \
        f'Expected regions file to have a contiguous set of indices from 0 to {largest_bin}.\n'
    
    try:
        subprocess.check_output(['sort', '-k','1,1', '-k','2,2n', '-c', regions_file])
    except subprocess.CalledProcessError:
        raise ValueError(
            f'Expected regions file to be sorted by chromosome and start position.\n'
            f'Use \'sort -k1,1 -k2,2n -o <filename> <filename>\' to sort the file.\n'
        )


def _make_fixed_size_windows(*, 
                            genome_file, 
                            window_size,
                            blacklist_file=None,
                            output = sys.stdout
                        ):
    
    process_kw = dict(
        universal_newlines=True,
        bufsize=10000,
    )

    makewindows_process = subprocess.Popen(
        ['bedtools', 'makewindows', '-g', genome_file, '-w', str(window_size)],
        stdout = subprocess.PIPE,
        **process_kw,
    )

    sort_process = subprocess.Popen(
        ['sort', '-k1,1', '-k2,2n'],
        stdin = makewindows_process.stdout,
        stdout = subprocess.PIPE,
        **process_kw,
    )

    if blacklist_file is not None:
        subract_process = subprocess.Popen(
            ['bedtools', 'intersect', '-a', '-', '-b', blacklist_file, '-v'],
            stdin = sort_process.stdout,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            bufsize=10000,
        )
        sort_process = subract_process

    add_id_process = subprocess.Popen(
        ['awk','-v','OFS=\t', '{print $0,NR-1,"0","+",$2,$3,"0,0,0","1",$3-$2,"0"}'],
        stdin = sort_process.stdout,
        stdout = output,
        **process_kw,
    )

    add_id_process.wait()


def stream_bedfile(bedfile):

    with open(bedfile, 'r') as f:

        for line in f:

            if line.startswith('#'):
                continue

            cols = line.strip().split('\t')

            if len(cols) < 3:
                raise ValueError(f'Bedfile {bedfile} must have at least 3 columns')
            
            feature = '1' if len(cols) == 3 else cols[3]

            chrom, start, end = cols[:3]
            start = int(start); end = int(end)

            yield chrom, start, end, feature


def _get_endpoints(*bedfiles):

    def _get_endpoints_bedfile(bedfile, track_id):
    
        for chrom, start, end, feature in stream_bedfile(bedfile):
            yield chrom, start, track_id, feature, True
            yield chrom, end, track_id, feature, False

    endpoints = [_get_endpoints_bedfile(bedfile, os.path.basename(bedfile)) for bedfile in bedfiles]
    
    return interleave_streams(*endpoints, key = lambda x : (x[0], x[1]))


def _endpoints_to_segments(endpoints): # change default min_windowsize 3 to 4

    active_features = Counter()
    feature_combination_ids = dict()
    prev_chrom = None; prev_pos = None

    for (chrom, pos, track_id, feature, is_start) in endpoints:

        pos = int(pos)

        if prev_chrom is None:
            prev_chrom = chrom; prev_pos = pos
        elif chrom != prev_chrom:
            active_features = Counter()
            prev_chrom = chrom; prev_pos = pos
        elif pos > prev_pos and len(active_features) > 0:
            
            is_nested_start = active_features[(track_id, feature)] > 0 and is_start
            is_nested_end = active_features[(track_id, feature)] > 1 and not is_start

            if is_nested_start or is_nested_end:
                '''logger.warning(
                    f'Nested features detected in file {track_id}, feature {feature} at position {chrom}:{pos}.\n'
                    'Avoided writing a new window.'
                )'''
                pass
            else:
                feature_combination = tuple(sorted(active_features.keys()))

                if not feature_combination in feature_combination_ids:
                    feature_combination_ids[feature_combination] = len(feature_combination_ids)

                yield chrom, prev_pos, pos, feature_combination_ids[feature_combination]

        if is_start:
            active_features[(track_id,feature)] += 1
        else:
            if active_features[(track_id,feature)] > 1:
                active_features[(track_id,feature)]-=1
            else:
                active_features.pop((track_id,feature))

        prev_pos = pos; prev_chrom = chrom


def format_bed12_record(region_id, segments):

    chrs, starts, ends = list(zip(*segments))
            
    region_start=min(starts); region_end=max(ends)
    region_chr = chrs[0]
    num_blocks = len(segments)
    
    block_sizes = ','.join(map(lambda x : str(x[0] - x[1]), zip(ends, starts)))
    block_starts = ','.join(map(lambda s : str(s - region_start), starts))

    return (
            region_chr,    # chr
            region_start,  # start
            region_end,    # end
            region_id,     # name
            '0','+',       # value, strand
            region_start,  # thickStart
            region_start,  # thickEnd
            '0,0,0',       # itemRgb,
            num_blocks,    # blockCount
            block_sizes,   # blockSizes
            block_starts,  # blockStarts
        )

def expand_args(func):
    def wrapper(x):
        return func(*x)
    return wrapper

def make_windows(
    *bedfiles,
    genome_file, 
    window_size, 
    blacklist_file, 
    output=sys.stdout, 
    min_windowsize=50,
):
    logger.info(f'Checking sort order of bedfiles ...')
    for bedfile in bedfiles:
        try:
            subprocess.run(["sort", "-k1,1", "-k2,2n", "--check", bedfile], check=True)
        except subprocess.CalledProcessError:
            raise ValueError(f'Bedfile {bedfile} is not sorted by chromosome and start position.')

    allowed_chroms=[]
    with open(genome_file,'r') as f:

        for line in f:
            if line.startswith('#'):
                continue
            allowed_chroms.append(line.strip().split('\t')[0].strip())

    logger.info(f'Using chromosomes: {", ".join(allowed_chroms)}')
    window_sizes=[]
    n_windows_written=[0]

    def accumulate_windowsizes(r):
        window_sizes.append(sum(map(lambda s : s[2] - s[1], r)))
        n_windows_written[0]+=1
        if n_windows_written[0] % 10000 == 0:
            logger.info(f'Wrote {n_windows_written[0]} windows ...')
        return r

    with tempfile.NamedTemporaryFile('w') as windows_file:

        logger.info(f'Making initial coarse-grained regions ...')
        _make_fixed_size_windows(
            genome_file=genome_file,
            window_size=window_size,
            output=windows_file,
        )
        windows_file.flush()
        
        '''
        Welllll..... what's going on here? I promise, it's beautiful in principle
        and only ugly in practice. If python had better syntax for functional programming
        and streaming iterators, this would look like a nice clean pipeline. Here it is in
        a more straightforward form:

        get_endpoints(*bedfiles) => map(filter_chroms(allowed_chroms)) => collect_to_segments \
            => map(convert_to_RegionOverlapComparitor) => filter_intersection(blacklist_file) \
            => map(unwrap_from_RegionOverlapComparitor) => collect_to_regions \
            => map(chop_to<chrom,start,end>) \
            => map(filter(min_windowsize)) => enumerate => map(format_bed12_record) \
            => write_to_output

        So why do it this way? Because this massive transformation is defined in terms of 
        lazy iterators, it's memory overhead is very small. This is important because 
        it means you can use an obscenely dense feature landscape without running out of memory.
        '''
        logger.info(f'Segmenting genome based on feature tracks ...')
        for bed12_record in \
                    map(expand_args(format_bed12_record),
                    enumerate(
                    map(accumulate_windowsizes, # this is here just to count the windows as they pass by. No transformation occurs.
                    filter(lambda r : sum(map(lambda s : s[2] - s[1], r)) > min_windowsize,
                    map(lambda r : list(map(lambda s : s[:3], r)),
                    buffered_aggregator(
                        map(RegionOverlapComparitor.unwrap,
                            filter_intersection(
                                map(expand_args(RegionOverlapComparitor),
                                    _endpoints_to_segments(
                                        filter(
                                            lambda x : x[0] in allowed_chroms,
                                            _get_endpoints(
                                                windows_file.name,
                                                *bedfiles,
                                            )
                                        )
                                    )
                                ),    
                                map(expand_args(RegionOverlapComparitor),
                                    stream_bedfile(blacklist_file)
                                ),
                            ),
                        ),
                        has_lapsed = lambda x, y : x[0] != y[0] or x[1] - y[1] > window_size*2,
                        key = lambda x : x[3]
                    )))))):
            print(*bed12_record, sep='\t', file=output)  

    q=(0.1, 0.25, 0.5, 0.75, 0.9)
    windowsize_dist=quantile(window_sizes, q)
    print(
f'''Window size report
------------------
Num windows   | {len(window_sizes)}
Smallest      | {min(window_sizes)}
Largest       | {max(window_sizes)}    
''' + '\n'.join(('Quantile={: <4} | {}'.format(str(k),str(int(v))) for k,v in zip(q, windowsize_dist)))
    )
