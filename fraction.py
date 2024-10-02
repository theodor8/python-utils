
def gcd(a, b):
    # a = q * b + r
    if b == 0: return a
    q = a // b
    r = a - q * b
    return gcd(b, r)



class Fraction:
    def __init__(self, num, den=1):
        if type(num) != int or type(den) != int:
            raise Exception("Not ints")
        d = gcd(num, den)
        self.num = int(num / d)
        self.den = int(den / d)

    @staticmethod
    def from_float(n):
        for i in range(1, 100000):
            if (n * i) % 1 == 0:
                return Fraction(n * i, i)

    def __str__(self):
        return str(self.num) if self.den == 1 else f'({self.num}/{self.den})'

    def eval(self):
        return self.num / self.den
    
    def reduced(self):
        return Fraction(self.num, self.den)
    
    def __add__(self, other):
        if type(other) == int:
            other = Fraction(other)
        return Fraction(self.num * other.den + other.num * self.den, self.den * other.den)
    
    def __radd__(self, other):
        return self + other
    
    def __sub__(self, other):
        return self + -other
    
    def __rsub__(self, other):
        return -self + other
    
    def __mul__(self, other):
        if type(other) == int:
            other = Fraction(other)
        return Fraction(self.num * other.num, self.den * other.den)

    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other):
        if type(other) == int:
            other = Fraction(other)
        return self * other.inv()
    
    def __rtruediv__(self, other):
        return other * self.inv()
    
    def __pow__(self, other):
        return Fraction(self.num ** other, self.den ** other)
    
    def __neg__(self):
        return Fraction(-self.num, self.den)
    
    def inv(self):
        return Fraction(self.den, self.num)
    
