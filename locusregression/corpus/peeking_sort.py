from collections import defaultdict

def sorted_iterator(iter, key = lambda x : x):
    # this is a generator that yields the next item in the iterator, 
    # but raises a ValueError if the items are not in order
    for curr in iter:

        try:
            _prev
        except UnboundLocalError:
            _prev = curr
        else:
            if not key(curr) >= key(_prev):
                raise ValueError('Items are out of order, {} should be greater than {}'\
                        .format(str(key(curr)), str(key(_prev)))
                    )
            
            _prev = curr

        yield curr


class PeekIterator:
    '''
    A class that wraps an iterator and allows you to peek at the next value,
    without advancing the iterator. It implements the Iterator trait
    but with the "peek" method added.
    '''
    def __init__(self, iterator):
        self._iterator = iterator
        self._depleted = False
        try:
            self._next = self._get_next()
        except StopIteration:
            self._depleted = True

    def _get_next(self):
        return next(self._iterator)
        
    def __next__(self):

        if self._depleted:
            raise StopIteration()            

        ret_value = self._next 
        
        try:
            self._next = self._get_next()
        except StopIteration:
            self._depleted = True
        
        return ret_value

    def peek(self):
        if self._depleted:
            raise StopIteration()

        return self._next

    def has_next(self):
        return not self._depleted

    def __eq__(self, other):
        return self.peek() == other.peek()

    def __gt__(self, other):
        return self.peek() > other.peek()
    

def interleave_streams(*iterators, key = lambda x : x):

    streams = [
        PeekIterator(sorted_iterator(stream, key=key))
        for stream in iterators
    ]

    while True:

        streams = [stream for stream in streams if stream.has_next()]
        
        if len(streams) == 0:
            break
        else:
            yield next(min(streams, key = lambda x : key(x.peek())))


def buffered_aggregator(
        iterator, 
        has_lapsed = lambda x, y : False, 
        key = lambda x : x,
        order_key = lambda x : x,
    ):
    '''
    This function takes an iterator and a function that determines whether a window has lapsed
    and returns an iterator that aggregates the values into windows. The function has_lapsed
    should take two arguments, the current value and the previous value, and return True if the
    window has lapsed.
    '''
    
    buffer = defaultdict(list)
    buffer_order = {}

    def get_lapsed_key(val):
        for dictkey, values in buffer.items():
            if has_lapsed(val, values[0]) \
                    and buffer_order[dictkey] == min(buffer_order.values()):
                return dictkey
        else:
            return None

    def pop_lapsed_keys(val):
        lapsed_key = get_lapsed_key(val)
        while not lapsed_key is None:
            yield lapsed_key
            lapsed_key = get_lapsed_key(val)

    while True:

        try:
            val = next(iterator)
        except StopIteration:
            break
        
        for dictkey in pop_lapsed_keys(val):
            buffer_order.pop(dictkey)
            yield buffer.pop(dictkey)
        
        dictkey=key(val)
        if not dictkey in buffer:
            buffer_order[dictkey] = order_key(val)
        buffer[dictkey].append(val)
        
    for dictkey in pop_lapsed_keys(val):
        buffer_order.pop(dictkey)
        yield buffer.pop(dictkey)


def filter_intersection(filter_iter, blacklist_iter):

    iter1 = PeekIterator(filter_iter)
    iter2 = PeekIterator(blacklist_iter)

    while iter1.has_next():

        while iter2.has_next() and iter1.peek() > iter2.peek():
            next(iter2)

        if iter2.has_next():
            if iter1.peek() == iter2.peek():
                next(iter1) # advance iter1, filter out the value
            else:
                yield next(iter1) # advance iter1 and yield the value
        else:
            yield next(iter1) # if there are no more values in iter2, 
                              # yield the value from iter1
