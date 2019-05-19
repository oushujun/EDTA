#!/usr/bin/env python
# -*- coding: utf-8 -*-


class cachedproperty(object):

    '''
    >>> class C(object):
    ...    def __init__(self, x):
    ...        self._x = x
    ...    def get_x(self):
    ...        print "getting x"
    ...        return self._x
    ...
    ...    def set_x(self, newx):
    ...        print "setting x with %s" % newx
    ...        self._x = newx
    ...    def del_x(self):
    ...        self._x = "i am deleted"
    ...
    ...    x = cachedproperty(get_x, set_x, del_x)
    ...    other_x = cachedproperty(get_x, set_x)

    >>> c = C(5)
    >>> c.x
    getting x
    5

    # cached.
    >>> c.x
    5

    >>> c.x, c.y = 6, 7
    setting x with 6

    # uncached.
    >>> c.x
    getting x
    6

    >>> c.x, c.y
    (6, 7)

    >>> c.y = 35
    >>> c.x, c.y
    (6, 35)

    # ok with multiple instances.
    >>> d = C(4)
    >>> d.x
    getting x
    4

    >>> c.x
    6
    >>> c.other_x = 7
    setting x with 7

    >>> c.get_x()
    getting x
    7
    >>> del c.x
    >>> c.x
    getting x
    \'i am deleted\'

    >>> c.set_x(22)
    setting x with 22

    # but the property cant know about it...
    >>> c.x
    \'i am deleted\'


    '''

    __slots__ = ("fget", "fset", "fdel", "n")

    def __init__(self, fget=None, fset=None, fdel=None):
        self.fget = fget
        self.fset = fset
        self.fdel = fdel
        self.n = "__" + fget.__name__

    def __get__(self, o, otype=None):
        if o is None:
            return None
        if self.n in o.__dict__:
            return (o.__dict__)[self.n]
        result = (o.__dict__)[self.n] = self.fget(o)
        return result

    def __set__(self, o, value):
        if self.fset is None:
            raise AttributeError("unsettable %s (with %s)" % (self.n, value))
        else:
            if self.n in o.__dict__:
                del (o.__dict__)[self.n]
            self.fset(o, value)

    def __delete__(self, o):
        if self.fdel is None:
            raise AttributeError("undeletable %s (with %s)" % (self.n, value))
        else:
            if self.n in o.__dict__:
                del (o.__dict__)[self.n]
            self.fdel(o)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
