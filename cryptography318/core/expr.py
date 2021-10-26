from numbers import *

from math import prod, gcd
from abc import ABCMeta, abstractmethod
from functools import reduce

from .tools import join_dict


class Expr(metaclass=ABCMeta):
    __slots__ = '_args'  # type: dict

    @property
    @abstractmethod
    def eval(self): ...

    @abstractmethod
    def simplify(self, *, _update=False): ...

    @property
    @abstractmethod
    def gcd(self): ...

    @abstractmethod
    def __str__(self): ...

    @abstractmethod
    def __repr__(self): ...

    @abstractmethod
    def __add__(self, other): ...

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    @abstractmethod
    def __mul__(self, other): ...

    def __rmul__(self, other):
        return self.__mul__(other)

    @abstractmethod
    def __truediv__(self, other): ...

    @abstractmethod
    def __rtruediv__(self, other): ...

    def __int__(self):
        return int(float(self))

    @abstractmethod
    def __float__(self): ...

    def __complex__(self):
        return complex(float(self))

    @abstractmethod
    def __abs__(self): ...

    @abstractmethod
    def __pos__(self): ...

    @abstractmethod
    def __neg__(self): ...

    @abstractmethod
    def __eq__(self, other): ...

    def __lt__(self, other):
        return self.eval < other

    def __le__(self, other):
        return self.eval <= other

    def __gt__(self, other):
        return self.eval > other

    def __ge__(self, other):
        return self.eval >= other

    @property
    def args(self: 'Expr'):
        return reduce(lambda r, c: r + tuple(c), self._args.values(), ())

    @property
    def dict(self: 'Expr'):
        return self._args

    @property
    def types(self: 'Expr'):
        return tuple(self._args.keys())

    @property
    def numerator(self):
        return +self

    @property
    def denominator(self):
        return 1


class Mul(Expr):

    def __init__(self: 'Mul', *args, _dict=None):

        # assumes construction is being performed inside Expr method, input is trusted
        if isinstance(_dict, dict):
            self._args = {key: value[:] for key, value in _dict.items()}
        else:
            self._args = {int: [1]}

        if args:
            for e in args:
                if isinstance(e, int) or (isinstance(e, float) and e.is_integer()):
                    self._args[int] = [prod(self._args.get(int, [1])) * int(e)]
                elif isinstance(e, Mul):
                    self._args = join_dict(self._args, e.dict, _null=[])
                elif type(e) in self._args:
                    self._args[type(e)].append(e)
                elif isinstance(e, Real):
                    self._args[type(e)] = [e]
                else:
                    raise TypeError(f"all arguments must inherit from Real")

        _int = prod(self._args.get(int, [1]))
        if not _int:
            self._args = {int: [0]}
        else:
            self._args[int] = [_int]

    @property
    def eval(self):
        return prod(float(e) for e in self.args)

    def simplify(self: 'Mul', *, _update=False):
        if len(self.args) < 2:
            return None if _update else +self
        else:
            reduced = {}
            # we know how to quickly reduce ints, so dont do that here
            _args = {key: value[:] for key, value in self._args.items() if key is not int}
            for t in _args:
                # print(f"reducing {t.__name__} ...")
                reduced = join_dict(reduce(lambda r, c: r * c, _args[t], Mul()).dict, reduced, _null=[])
                # print(f"{t.__name__} reduced: {({key.__name__: value for key, value in reduced.items()})}\n")

            # reduce the ints (including any new ints created from result of any Real products)
            reduced[int] = [prod(reduced.get(int, [1])) * prod(self._args.get(int, [1]))]

            # print(f"all reduced: {({key.__name__: value for key, value in reduced.items()})}")

            if _update:
                self._args = reduced
            else:
                return Mul(_dict=reduced)

    @property
    def gcd(self: 'Mul'):
        return prod(self._args.get(int, [1]))

    def print_dict(self: 'Mul'):
        print({key.__name__: value for key, value in self._args.items()})

    @staticmethod
    def _validate_div(_div, _args, _type_other, _type_num, index):
        """
        :param _div: result of division
        :param _args: dictionary of expression
        :param _type_other: type of divisor (not guaranteed to be a type in _args)
        :param _type_num: type of numerator in division (guaranteed to be a type in _args)
        :param index: index in list of _args[_type_num] of numerator
        """
        # most common if _div is int
        if isinstance(_div, float):
            # if float and not integer, definitely won't be a valid value later
            if _div.is_integer():
                del _args[_type_num][index]
                _args[int] = [prod(_args.get(int, [1])) * int(_div)]
                return Mul(_dict={key: value for key, value in _args.items() if value})
            else:
                return False
        # ints will pass with this in addition to Rational objects
        elif isinstance(_div, Rational) and _div.denominator == 1 and isinstance(_div.numerator, int):
            del _args[_type_num][index]
            _args[int] = [prod(_args.get(int, [1])) * _div.numerator]
            return Mul(_dict={key: value for key, value in _args.items() if value})
        # if same type as current type in iteration, replace a with new value
        elif isinstance(_div, _type_num):
            _args[_type_num][index] = _div
            return Mul(_dict=_args)
        # if same type as other, this is considered a success
        elif isinstance(_div, _type_other):
            # remove current instance a
            del _args[_type_num][index]
            # add result _div either to current list or empty list
            _args[_type_other] = _args.get(_type_other, [])
            _args[_type_other].append(_div)
            return Mul(_dict=_args)
        elif isinstance(_div, Mul):
            del _args[_type_num][index]
            _args = join_dict(_args, _div.dict, _null=[])
            return Mul(_dict={key: value for key, value in _args.items() if value})
        else:
            return False

    def __str__(self):
        if int in self._args:
            k = prod(self._args[int])
            self._args[int] = [k]
            if k == 1:
                del self._args[int]
                return '*'.join(str(e) for e in self.args)
            elif len(self.args) == 1:
                return str(k)
            else:
                return f'{k}*' + '*'.join(str(e) for e in self.args if not isinstance(e, int))
        else:
            return '*'.join(str(e) for e in self.args)

    def __repr__(self):
        return f"{self.__class__.__name__}{self.args}"

    def __add__(self: 'Mul', other):
        if isinstance(other, Mul):
            # if equal terms except for the ints then add the ints and return same term
            # ex. 2*Sqrt(3)*Sqrt(2) + 5*Sqrt(3)*Sqrt(2) = 8*Sqrt(3)*Sqrt(2)
            if sorted(e for e in other.args if not isinstance(e, int)) == \
                    sorted(e for e in self.args if not isinstance(e, int)):
                _args = {key: value for key, value in self._args.items()}
                _args[int] = [prod(_args.get(int, [1])) + prod(other.dict.get(int, [1]))]
                return Mul(_dict=_args)
            else:
                return Add(self, other)
        elif isinstance(other, Add):
            _args = {key: value[:] for key, value in other.dict.items()}
            if Mul in other.dict:
                for i, e in enumerate(_args[Mul]):
                    _try = self + e
                    if isinstance(_try, Mul):
                        _args[Mul][i] = _try
                        return Add(_dict=_args)
                _args[Mul].append(self)
                return Add(_dict=_args)
            else:
                _args[Mul] = [self]
                return Add(_dict=_args)
        elif isinstance(other, Real) and not isinstance(other, float):
            _args = {key: value for key, value in self._args.items() if key is not int}
            # if only one type of arg; if only one arg in that type;
            if len(_args) == 1 and len(v := tuple(_args.values())[0]) == 1:
                # if the only arg in the Mul directly equals other, then increment the integer value of Mul and return
                if v[0] == other:
                    _args[int] = [prod(self._args.get(int, [1])) + 1]
                    return Mul(_dict=_args)
            # anything other than that, create an Add expr
            return Add(self, other)
        else:
            return NotImplemented

    def __mul__(self: 'Mul', other):
        # if Mul(0), should immediately return Mul(0)
        if not (_int := prod(self._args.get(int, [1]))):
            return Mul(_dict={int: [0]})
        elif isinstance(other, int):
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy
            _args[int] = [prod(_args.get(int, [1])) * other]
            return Mul(_dict=_args)
        elif isinstance(other, Mul):
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy

            for _type in other.dict:
                if _type is int:
                    _args[int] = [prod(_args.get(int, [1])) * prod(other.dict[int])]
                elif _type in _args:
                    for a in other.dict[_type]:
                        added = False
                        vals = _args[(t := type(a))]
                        for i, v in enumerate(vals):
                            _try = a * v
                            if isinstance(_try, Rational) and _try.denominator == 1 and isinstance(_try.numerator, int):
                                _args[int] = [prod(_args.get(int, [1])) * _try.numerator]
                                del _args[t][i]
                                added = True
                                break
                            elif isinstance(_try, int):
                                _args[int] = [prod(_args.get(int, [1])) * _try]
                                del _args[t][i]
                                added = True
                                break
                            elif isinstance(_try, t):
                                _args[t][i] = _try
                                added = True
                                break
                        # if unable to take product with any other elements of same type, add it to list for that type
                        if not added:
                            _args[t].append(a)
                else:
                    _args[_type] = other.dict[_type]

            # if value is empty list, then all of that type was simplified to int, so remove it
            return Mul(_dict={key: value for key, value in _args.items() if value})
        elif type(other) in self._args:
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy
            added = False
            vals = _args[(t := type(other))]
            for i, v in enumerate(vals):
                _try = other * v
                if isinstance(_try, Rational) and _try.denominator == 1 and isinstance(_try.numerator, int):
                    _args[int] = [prod(_args.get(int, [1])) * _try.numerator]

                    if len(_args[t]) == 1:
                        del _args[t]
                    else:
                        del _args[t][i]
                    added = True
                    break
                elif isinstance(_try, int):
                    _args[int] = [prod(_args.get(int, [1])) * _try]

                    if len(_args[t]) == 1:
                        del _args[t]
                    else:
                        del _args[t][i]
                    added = True
                    break
                elif isinstance(_try, t):
                    _args[t][i] = _try
                    added = True
                    break
                # if Mul, remove original instance then add all into new Mul expr
                elif isinstance(_try, Mul):
                    if len(_args[t]) == 1:
                        del _args[t]
                    else:
                        del _args[t][i]
                    return Mul(*reduce(lambda r, c: r + tuple(c), _args.values(), ()), _try)

            if not added:
                _args[t].append(other)

            return Mul(_dict=_args)
        elif isinstance(other, Real) and not isinstance(other, float):
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy
            _args[type(other)] = [other]

            return Mul(_dict=_args)
        else:
            return NotImplemented

    def __truediv__(self: 'Mul', other):
        if isinstance(other, int):
            _args = {key: value[:] for key, value in self._args.items()}

            # create _types s.t. if int is an arg, it will come first since that is most likely to be a value
            # that produces a valid result, after that, continue to loop through the rest of the types
            if int in _args:
                _types = [int] + [t for t in self.types if t is not int]
                _args[int] = [prod(_args[int])]
            else:
                _types = [*self.types]

            for _t in _types:
                for i, e in enumerate(_args[_t]):
                    try:
                        _div = e / other
                    except ValueError or TypeError:
                        raise ValueError(f"{self} is not divisible by {other}")
                    else:
                        res = self._validate_div(_div, _args, int, _t, i)
                        if res:
                            return res
            # if never returned in loop, nothing produces a valid result so this value does not divide the expression
            raise ValueError(f"{self} is not divisible by {other}")
        elif isinstance(other, Mul):
            # list of args without int values
            no_ints = lambda d: filter(lambda e: not isinstance(e, int), d.args)
            # if args of both Mul's are the same and integer values of other divide integers of self, return new Mul
            if sorted(no_ints(self)) == sorted(no_ints(other)) and \
                    (s := prod(self._args.get(int, [1]))) % (o := prod(other.dict.get(int, [1]))) == 0:
                return Mul(s // o)
            # if they aren't the same terms (exluding ints) then just try division at each step, returning
            # result if no errors are thrown
            else:
                try:
                    _div = +self
                    for e in other.args:
                        _div /= e
                except ValueError:
                    raise ValueError(f"{self} is not divisible by {other}")
                else:
                    return _div
        elif isinstance(other, Real) and not isinstance(other, float):
            _args = {key: value[:] for key, value in self._args.items()}
            if other in self._args.get((t := type(other)), []):
                _args[t].remove(other)
                return Mul(_dict=_args)
            else:
                # create list of types so that if type(other) is in self._args, it comes first since it is
                # most likely to produce a valid output, then iterate through the rest of the types
                if t in self._args:
                    _types = [t] + [e for e in self.types if e is not t]
                else:
                    _types = [e for e in self._args.keys()]

                # the goal is to try every combination of division, if any result in an output we like, use it
                # if we don't get any valid outputs, then we cannot divide
                for _t in _types:
                    for i, a in enumerate(_args[_t]):
                        try:
                            _div = a / other
                        except TypeError or ValueError:
                            continue
                        else:
                            res = self._validate_div(_div, _args, t, _t, i)
                            # res is either False or a Mul instance
                            if res:
                                return res
                raise ValueError(f"{self} is not divisible by {other}")
        elif isinstance(other, Add):
            if len(other.args) == 1 and isinstance(other.args[0], Mul):
                return self / other.args[0]
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, (int, Real, Add)) and not isinstance(other, float):
            _div = other
            try:
                for e in self.args:
                    _div /= e
            except TypeError or ValueError:
                raise ValueError(f"{other} is not divisible by {self}")
            else:
                return Mul(_div)
        else:
            return NotImplemented

    def __float__(self):
        return prod(float(e) for e in self.args)

    def __bool__(self: 'Mul'):
        # if self = Mul(0) return False right away, otherwise check self.eval
        if prod(self._args[int]) == 0 or not self.eval:
            return False
        else:
            return True

    def __hash__(self):
        return hash(self.eval)

    def __abs__(self):
        if self.eval < 0:
            return -self
        else:
            return +self

    def __pos__(self: 'Mul'):
        return Mul(_dict=self._args)

    def __neg__(self: 'Mul'):
        _args = {key: value[:] for key, value in self._args.items()}
        if int in self._args:
            _args[int] = [-prod(_args[int])]
            return Mul(_dict=_args)
        else:
            _args[int] = [-1]
            return Mul(_dict=_args)

    def __eq__(self, other):
        if isinstance(other, Mul):
            return sorted(self.args) == sorted(other.args)
        else:
            return self.eval == other


class Add(Expr):

    def __init__(self: 'Add', *args, _dict=None):
        if isinstance(_dict, dict):
            self._args = {key: value[:] for key, value in _dict.items()}
        else:
            self._args = {int: [0]}

        if args:
            for e in args:
                if isinstance(e, int) or (isinstance(e, float) and e.is_integer()):
                    self._args[int] = [sum(self._args.get(int, [0])) + int(e)]
                elif isinstance(e, Add):
                    self._args = join_dict(self._args, e.dict, _null=[])
                elif isinstance(e, Mul) and len(e.dict) == 1 and int in e.dict:
                    self._args[int] = [sum(self._args.get(int, [0])) + prod(e.dict[int])]
                elif type(e) in self._args:
                    self._args[type(e)].append(e)
                elif isinstance(e, Mul):
                    self._args[Mul] = [e]
                elif isinstance(e, Real) and not isinstance(e, float):
                    self._args[type(e)] = [e]
                else:
                    raise TypeError(f"all arguments must inherit from Real")

    @property
    def eval(self):
        return sum(float(e) for e in self.args)

    @property
    def args(self: 'Add'):
        # integer value of 0 should not be returned as an arg (prevents ZeroDivision)
        if int in self._args:
            k = sum(self._args[int])
            if len(self._args) == 1:
                return tuple([k])
            elif k == 0:
                del self._args[int]
            else:
                self._args[int] = [k]
        return reduce(lambda r, c: r + tuple(c), self._args.values(), ())

    def simplify(self: 'Add', *, _update=False):
        # remove any null values regardless of _update's status
        self._args = {key: list(e for e in value if e) for key, value in self._args.items()}
        if len(self.args) < 2:
            return None if _update else +self
        else:
            reduced = {}
            _args = {key: value[:] for key, value in self._args.items()}
            for t in _args:
                reduced = join_dict(reduce(lambda r, c: r + c, _args[t], Add()).dict, reduced, _null=[])

            # because some types can multiply to result in ints, join dict might return reduced w/ multiple
            # entries in int, so take product before returning
            reduced[int] = [sum(reduced.get(int, [0]))]

            if _update:
                self._args = reduced
            else:
                return Add(_dict=reduced)

    @property
    def gcd(self: 'Add'):
        if int in self._args:
            if len(self._args) == 1:
                k = sum(self._args[int])
                self._args[int] = [k]
                return k
            elif len(self._args) == 2 and Mul in self._args:
                k = sum(self._args[int])
                g = gcd(k, *map(lambda e: e.gcd, self._args[Mul]))
                self._args[int] = [k]
                return g
            else:
                return 1
        elif len(self._args) == 1 and Mul in self._args:
            return gcd(*map(lambda e: e.gcd, self._args[Mul]))
        else:
            return 1

    def print_dict(self: 'Add'):
        print({key.__name__: value for key, value in self._args.items()})

    def __str__(self):
        if int in self._args:
            k = sum(self._args[int])
            self._args[int] = [k]
            if k == 0:
                del self._args[int]
                return '+'.join(str(e) for e in self.args)
            elif len(self.args) == 1:
                return str(k)
            else:
                return f'{k}+' + '+'.join(str(e) for e in self.args if not isinstance(e, int))
        else:
            return '+'.join(str(e) for e in self.args)

    def __repr__(self):
        return f"{self.__class__.__name__}{self.args}"

    def __add__(self: 'Add', other):
        if isinstance(other, int):
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy
            _args[int] = [sum(_args.get(int, [0])) + other]
            return Add(_dict=_args)
        elif isinstance(other, Add):
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy

            for _type in other.dict:
                if _type is int:
                    _args[int] = [sum(_args.get(int, [0])) + sum(other.dict[int])]
                elif _type in _args:
                    for a in other.dict[_type]:
                        added = False
                        vals = _args[(t := type(a))]
                        for i, v in enumerate(vals):
                            if a == v:
                                if Mul in _args:
                                    _args[Mul].append(Mul(2, a))
                                else:
                                    _args[Mul] = [Mul(2, a)]
                                del _args[t][i]
                                added = True
                                break
                            _try = a + v
                            if isinstance(_try, Rational) and _try.denominator == 1 and isinstance(_try.numerator, int):
                                _args[int] = [sum(_args.get(int, [0])) + _try.numerator]
                                del _args[t][i]
                                added = True
                                break
                            elif isinstance(_try, int):
                                _args[int] = [sum(_args.get(int, [0])) + _try]
                                del _args[t][i]
                                added = True
                                break
                            elif isinstance(_try, (t, Mul)):
                                _args[t][i] = _try
                                added = True
                                break

                        # if unable to take product with any other elements of same type, add it to list for that type
                        if not added:
                            _args[t].append(a)
                else:
                    _args[_type] = other.dict[_type]

            # if value is empty list, then all of that type was simplified to int, so remove it
            return Add(_dict={key: value for key, value in _args.items() if value})
        elif type(other) in self._args:
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy
            added = False
            vals = _args[(t := type(other))]
            for i, v in enumerate(vals):
                if other == v:
                    if Mul in _args:
                        _args[Mul].append(Mul(2, other))
                    else:
                        _args[Mul] = [Mul(2, other)]
                    del _args[t][i]
                    added = True
                    break
                _try = other + v
                if isinstance(_try, Rational) and _try.denominator == 1 and isinstance(_try.numerator, int):
                    _args[int] = [sum(_args.get(int, [0])) + _try.numerator]
                    del _args[t][i]
                    added = True
                    break
                elif isinstance(_try, int):
                    _args[int] = [sum(_args.get(int, [0])) + _try]
                    del _args[t][i]
                    added = True
                    break
                elif isinstance(_try, (t, Mul)):
                    _args[t][i] = _try
                    added = True
                    break

            # if unable to take product with any other elements of same type, add it to list for that type
            if not added:
                _args[t].append(other)

            return Add(_dict={key: value for key, value in _args.items() if value})
        elif isinstance(other, Real) and not isinstance(other, float):
            _args = {key: value[:] for key, value in self._args.items()}  # dict copy
            _args[type(other)] = [other]
            return Add(_dict=_args)
        elif isinstance(other, Mul):
            # addition is commutative so just use same method that Mul has
            return other.__add__(self)
        else:
            return NotImplemented

    def __mul__(self: 'Add', other):
        self.simplify(_update=True)
        if isinstance(other, int):
            # skips a lot of attempted multiplication
            if other == 1:
                return +self
            # also skips a lot of multiplication, by using each Real's __neg__() method
            elif other == -1:
                return -self
            else:
                _args = {key: list(e * other for e in value) for key, value in self._args.items()}
                return Add(_dict=_args)
        elif isinstance(other, Add):
            return Add(*map(lambda e: self * e, other.args))
        elif isinstance(other, Real) and not isinstance(other, float):
            _try = [e * other for e in self.args]

            # maps indices, returning result of direct multiplication if a valid type, otherwise
            # creating a Mul object
            _args = map(lambda v: _try[v] if type(_try[v]) in (*self.types, Mul) else Mul(
                self.args[v], other), range(len(self.args)))
            return Add(*_args).simplify()
        elif isinstance(other, Mul):
            return reduce(lambda r, c: r * c, other.args, self)
        else:
            return NotImplemented

    def __truediv__(self: 'Add', other):
        if isinstance(other, int):
            if other == 1:
                return +self
            elif other == -1:
                return -self
            elif other == 0:
                raise ZeroDivisionError
            elif len(self._args) == 1 and int in self._args:
                d, r = divmod(sum(self._args[int]), other)
                if r == 0:
                    return Add(d)
                else:
                    raise ValueError(f"{self} is not divisible by {other}")
            else:
                try:
                    _args = []
                    for e in self.args:
                        _args.append(e / other)
                except TypeError or ValueError:
                    raise ValueError(f"{self} is not divisible by {other}")
                else:
                    return Add(*_args)
        elif isinstance(other, Mul):
            # if its a Mul, it can only divide an Add if it only has one expression
            try:
                if len(self.args) == 1:
                    return Add(self.args[0] / other)
                else:
                    raise ValueError
            except ValueError or TypeError:
                raise ValueError(f"{self} is not divisible by {other}")

        elif isinstance(other, Add):
            # two ways an Add expression can divide another: either they are equal, in which case the divison = 1
            # or the numerator expression is a multiple of the denominator, meaning when iterating through
            # sorted(args) of both, every division should return the same value
            if sorted(self.args) == sorted(other.args):
                return Add(1)
            elif len(self.args) == len(other.args):
                try:
                    _self, _other = sorted(self.args), sorted(other.args)
                    _int = _self[0] / _other[0]

                    if len(_self) == 1:
                        return Add(_int)
                    else:
                        for m1, m2 in zip(_self[1:], _other[1:]):
                            res = m1 / m2
                            if res != _int:
                                raise ValueError
                except ValueError or TypeError:
                    raise ValueError(f"{self} is not divisible by {other}")
                else:
                    return Add(_int)
            else:
                raise ValueError(f"{self} is not divisible by {other}")
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, (int, Mul, Real)) and not isinstance(other, float):
            try:
                if len(self.args) == 1:
                    return Add(other / self.args[0])
                else:
                    raise ValueError
            except ValueError or TypeError:
                raise ValueError(f"{other} is not divisible by {self}")
        else:
            return NotImplemented

    def __float__(self):
        return sum(float(e) for e in self.args)

    def __bool__(self):
        return bool(self.eval)

    def __abs__(self):
        if self.eval < 0:
            return -self
        else:
            return +self

    def __pos__(self: 'Add'):
        # Add.__init__() makes copies of value lists already
        return Add(_dict=self._args)

    def __neg__(self: 'Add'):
        _args = {key: list(-e for e in value) for key, value in self._args.items()}
        return Add(_dict=_args)

    def __eq__(self, other):
        if isinstance(other, Add):
            return sorted(self.args) == sorted(other.args)
        else:
            return self.eval == other
