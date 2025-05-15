from abc import ABC, abstractmethod
from math import comb, erf, exp, factorial, pi, sqrt, gamma
from calc import Func, Pol, Const, Var


class Distribution(ABC):

    @property
    @abstractmethod
    def E(self) -> float: # expected
        pass

    @property
    @abstractmethod
    def V(self) -> float: # variance
        pass

    @property
    def D(self) -> float:
        return sqrt(self.V)

    @abstractmethod
    def f(self, x: float) -> float: # mass/density function (PMF/PDF)
        pass

    @abstractmethod
    def F(self, x: float) -> float: # cumulative distribution function (CDF)
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass



class Disc(Distribution): # discrete
    def __init__(self, p: dict[float, float]):
        self.p = p

    @property
    def E(self) -> float:
        return sum(x * self.p[x] for x in self.p)

    @property
    def V(self) -> float:
        return Disc({x * x: self.p[x] for x in self.p}).E - self.E**2

    def f(self, x: float) -> float:
        if x not in self.p:
            raise ValueError
        return self.p[x]

    def F(self, x: float) -> float:
        return sum(self.p[xi] for xi in self.p if xi <= x)

    def __repr__(self):
        return f"Disc({self.p})"

class Cont(Distribution): # continuous
    def __init__(self, f: Func, a: float = float('-inf'), b: float = float('inf')):
        self.func = f
        self.a = a
        self.b = b

    @property
    def E(self) -> float:
        return (Var() * self.func).integral(self.a, self.b)

    @property
    def V(self) -> float:
        return Cont((Var() * self.func), self.a, self.b).E - self.E**2

    def f(self, x: float) -> float:
        if x < self.a or x > self.b:
            return 0
        return self.func.eval(x)

    def F(self, x: float) -> float:
        if x < self.a:
            return 0
        if x > self.b:
            return 1
        return self.func.integral(self.a, x)

    def __repr__(self):
        return f"Cont({self.func}, {self.a}, {self.b})"



class Bin(Distribution):
    def __init__(self, n: int, p: float):
        if n < 0 or p <= 0 or p > 1:
            raise ValueError
        self.n = n
        self.p = p

    @property
    def E(self) -> float:
        return self.n * self.p

    @property
    def V(self) -> float:
        return self.n * self.p * (1 - self.p)

    def f(self, x: float) -> float:
        if x % 1 != 0 or x < 0 or x > self.n:
            raise ValueError
        x = int(x)
        return comb(self.n, x) * self.p**x * (1 - self.p)**(self.n - x)

    def F(self, x: float) -> float:
        if x % 1 != 0 or x < 0 or x > self.n:
            raise ValueError
        x = int(x)
        return sum(self.f(xi) for xi in range(x + 1))
    
    def __add__(self, other):
        if not isinstance(other, Bin):
            raise TypeError
        if self.p != other.p:
            raise ValueError
        return Bin(self.n + other.n, self.p)

    def __repr__(self):
        return f"Bin({self.n}, {self.p})"


class Po(Distribution):
    def __init__(self, e: float):
        if e <= 0:
            raise ValueError
        self.e = e

    @property
    def E(self) -> float:
        return self.e

    @property
    def V(self) -> float:
        return self.e

    def f(self, x: float) -> float:
        if x % 1 != 0 or x < 0:
            raise ValueError
        x = int(x)
        return self.e**x * exp(-self.e) / factorial(x)

    def F(self, x: float) -> float:
        if x % 1 != 0 or x < 0:
            raise ValueError
        k = int(x)
        return exp(-self.e) * sum(self.e**j / factorial(j) for j in range(k + 1))
    
    def __add__(self, other):
        if not isinstance(other, Po):
            raise TypeError
        return Po(self.e + other.e)

    def __repr__(self):
        return f"Po({self.e})"

class Re(Distribution):
    def __init__(self, a: float, b: float):
        if a >= b:
            raise ValueError
        self.a = a
        self.b = b

    @property
    def E(self) -> float:
        return (self.a + self.b) / 2

    @property
    def V(self) -> float:
        return (self.b - self.a)**2 / 12

    def f(self, x: float) -> float:
        if x < self.a or x > self.b:
            return 0
        return 1 / (self.b - self.a)

    def F(self, x: float) -> float:
        if x < self.a:
            return 0
        if x > self.b:
            return 1
        return (x - self.a) / (self.b - self.a)

    def __repr__(self):
        return f"Re({self.a}, {self.b})"

class Exp(Distribution):
    def __init__(self, a: float):
        if a <= 0:
            raise ValueError
        self.a = a

    @property
    def E(self) -> float:
        return self.a

    @property
    def V(self) -> float:
        return self.a**2

    def f(self, x: float) -> float:
        if x < 0:
            raise ValueError
        return 1 / self.a * exp(-x / self.a)

    def F(self, x: float) -> float:
        if x < 0:
            return 0
        return 1 - exp(-x / self.a)

    def __repr__(self):
        return f"Exp({self.a})"
    

class N(Distribution):
    def __init__(self, e: float, v: float):
        if v < 0:
            raise ValueError
        self.e = e
        self.v = v

    @property
    def E(self) -> float:
        return self.e

    @property
    def V(self) -> float:
        return self.v

    def f(self, x: float) -> float:
        return 1 / (sqrt(2*pi*self.v)) * exp(-(x-self.e)**2/(2*self.v**2))

    def F(self, x: float) -> float:
        return 0.5 * (1 + erf((x - self.e) / sqrt(2 * self.v)))

    def __add__(self, other):
        if not isinstance(other, N):
            raise TypeError
        return N(self.e + other.e, self.v + other.v)

    def __sub__(self, other):
        if not isinstance(other, N):
            raise TypeError
        return N(self.e - other.e, self.v + other.v)

    def __repr__(self):
        return f"N({self.e}, {self.v})"





class Sample:
    def __init__(self, *x: float):
        self.x = x

    @property
    def E(self) -> float:
        return sum(self.x) / len(self.x)

    @property
    def V(self) -> float:
        n = len(self.x)
        xm = self.E
        return 1 / (n - 1) * sum((xi - xm)**2 for xi in self.x)

    @staticmethod
    def S(x: 'Sample', y: 'Sample') -> float:
        return sum((xi - x.E) * (yi - y.E) for xi, yi in zip(x.x, y.x))

    def __repr__(self):
        return f"Sample{self.x}"

    def __getitem__(self, i: int) -> float:
        return self.x[i]

    def __len__(self) -> int:
        return len(self.x)





class SamplePairs:
    def __init__(self, *x):
        if len(x) == 2 and isinstance(x[0], Sample) and isinstance(x[1], Sample):
            if len(x[0]) != len(x[1]):
                raise ValueError("Samples must be of the same length")
            self.x = x[0]
            self.y = x[1]
            return
        self.x = Sample(*[x[0] for x in x])
        self.y = Sample(*[x[1] for x in x])

    def __repr__(self):
        pairs = [f"({self.x[i]}, {self.y[i]})" for i in range(len(self.x))]
        return f"SamplePairs({', '.join(pairs)})"

    @property
    def lm(self) -> tuple[float, float]: # linear model, kx + m
        k = Sample.S(self.x, self.y) / Sample.S(self.x, self.x)
        m = self.y.E - k * self.x.E
        return k, m

    @property
    def r(self) -> float: # pearson correlation coefficient
        return Sample.S(self.x, self.y) / sqrt(Sample.S(self.x, self.x) * Sample.S(self.y, self.y))

    @property
    def s(self) -> float:
        return sqrt(1/(len(self.x) - 2) * (Sample.S(self.y, self.y) - Sample.S(self.x, self.y)**2 / Sample.S(self.x, self.x)))


