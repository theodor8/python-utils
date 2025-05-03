from abc import ABC, abstractmethod
from math import exp



class Func(ABC):

    def __add__(self, other: 'Func') -> 'Func':
        return Add(self, other)

    def __sub__(self, other: 'Func') -> 'Func':
        return Sub(self, other)

    def __mul__(self, other: 'Func') -> 'Func':
        return Mul(self, other)

    def __truediv__(self, other: 'Func') -> 'Func':
        return Div(self, other)

    def __pow__(self, other: 'Func') -> 'Func':
        return Pow(self, other)

    @abstractmethod
    def eval(self, x: float) -> float:
        pass

    @abstractmethod
    def der(self) -> 'Func': # derivative
        pass

    @abstractmethod
    def prim(self) -> 'Func': # antiderivative, primitive function
        pass

    def integral(self, a: float, b: float) -> float:
        p = self.prim()
        return p.eval(b) - p.eval(a)



class Unary(Func):
    def __init__(self, f: Func) -> None:
        self.f = f

class Binary(Func):
    def __init__(self, f1: Func, f2: Func) -> None:
        self.f1 = f1
        self.f2 = f2



class Add(Binary):
    def eval(self, x: float) -> float:
        return self.f1.eval(x) + self.f2.eval(x)

    def der(self) -> Func:
        return Add(self.f1.der(), self.f2.der())

    def prim(self) -> Func:
        return Add(self.f1.prim(), self.f2.prim())

class Sub(Binary):
    def eval(self, x: float) -> float:
        return self.f1.eval(x) - self.f2.eval(x)

    def der(self) -> Func:
        return Sub(self.f1.der(), self.f2.der())

    def prim(self) -> Func:
        return Sub(self.f1.prim(), self.f2.prim())

class Mul(Binary):
    def eval(self, x: float) -> float:
        return self.f1.eval(x) * self.f2.eval(x)

    def der(self) -> Func:
        return Add(Mul(self.f1.der(), self.f2), Mul(self.f1, self.f2.der()))

    def prim(self) -> Func:
        raise NotImplementedError

class Div(Binary):
    def eval(self, x: float) -> float:
        return self.f1.eval(x) / self.f2.eval(x)

    def der(self) -> Func:
        return Div(Sub(Mul(self.f1.der(), self.f2), Mul(self.f1, self.f2.der())), Mul(self.f2, self.f2))

    def prim(self) -> Func:
        raise NotImplementedError

class Pow(Binary):
    def eval(self, x: float) -> float:
        return self.f1.eval(x) ** self.f2.eval(x)

    def der(self) -> Func:
        raise NotImplementedError

    def prim(self) -> Func:
        raise NotImplementedError



class Exp(Unary):
    def eval(self, x: float) -> float:
        return exp(self.f.eval(x))

    def der(self) -> Func:
        return Mul(Exp(self.f), self.f.der())

    def prim(self) -> Func:
        raise NotImplementedError

class Pol(Func):
    def __init__(self, pol: dict[float, float]) -> None:
        self.pol = pol # key = exponent, value = coefficient

    def coef(self, exp: int) -> float: # if not sure if exp exists in dict
        if exp in self.pol: return self.pol[exp]
        return 0
    
    def is_null(self) -> bool:
        return all(c == 0 for c in self.pol.values())

    def der(self) -> 'Pol':
        result = {}
        for e, c in self.pol.items():
            if e == 0: continue
            result[e - 1] = c * e
        return Pol(result)

    def prim(self) -> 'Pol':
        result = {}
        for e, c in self.pol.items():
            result[e + 1] = c / (e + 1)
        return Pol(result)
    
    def eval(self, x: float) -> float:
        return sum(self.pol[e] * x**e for e in self.pol)

    def __str__(self) -> str:
        result = ''
        s = dict(sorted(self.pol.items(), key=lambda item: item[0], reverse=True))
        for e, c in s.items():
            if c == 0: continue
            if len(result) > 0 or c < 0: result += ' + ' if c > 0 else ' - '
            if abs(c) != 1 or e == 0: result += str(abs(c))
            if e != 0: result += 'x' if e == 1 else f'x^{e}'
        return result

    def __mul__(self, other):
        if not isinstance(other, Pol):
            return super().__mul__(other)
        result = {}
        for e1, c1 in self.pol.items():
            for e2, c2 in other.pol.items():
                e = e1 + e2
                c = c1 * c2
                if e not in result: result[e] = 0
                result[e] += c
        return Pol(result)

    def __add__(self, other):
        if not isinstance(other, Pol):
            return super().__mul__(other)
        result = self.pol.copy()
        for e, c in other.pol.items():
            if e not in result: result[e] = 0
            result[e] += c
        return Pol(result)

    def __sub__(self, other):
        if not isinstance(other, Pol):
            return super().__mul__(other)
        result = self.pol.copy()
        for e, c in other.pol.items():
            if e not in result: result[e] = 0
            result[e] -= c
        return Pol(result)

    def __pow__(self, other):
        if not isinstance(other, int):
            return super().__mul__(other)
        result = Pol({0: 1})
        for _ in range(other):
            result *= self
        return result


class Const(Pol):
    def __init__(self, c: float) -> None:
        super().__init__({0: c})


class Var(Pol):
    def __init__(self) -> None:
        super().__init__({1: 1})

