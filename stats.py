from abc import ABC, abstractmethod
from math import comb


class Distribution(ABC):

    @property
    @abstractmethod
    def E(self): # expected
        pass

    @property
    @abstractmethod
    def V(self): # variance
        pass

    @property
    @abstractmethod
    def f(self): # mass/density function (PMF/PDF)
        pass

    @property
    @abstractmethod
    def F(self): # cumulative distribution function (CDF)
        pass


class Bin(Distribution):
    def __init__(self, n, p):
        if n < 0 or p <= 0 or p > 1:
            raise ValueError
        self.n = n
        self.p = p

    @property
    def E(self):
        return self.n * self.p

    @property
    def V(self):
        return self.n * self.p * (1 - self.p)

    @property
    def f(self):
        def f(k):
            if k % 1 != 0 or k < 0 or k > self.n:
                raise ValueError
            return comb(self.n, k) * self.p**k * (1 - self.p)**(self.n - k)
        return f

    @property
    def F(self):
        def F(k):
            if k % 1 != 0 or k < 0 or k > self.n:
                raise ValueError
            return sum(self.f(ki) for ki in range(k + 1))
        return F


a = Bin(2, 0.5)
print(a.f(0.0))

