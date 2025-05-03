
class Function:
    def __init__(self, f):
        self.f = f
        self.dx = 0.00001
        self.a = -10
        self.b = 10

    def value(self, x):
        return self.f(x)

    def __add__(self, other):
        return Function(lambda x: self.value(x) + other.value(x))
    
    def __sub__(self, other):
        return Function(lambda x: self.value(x) - other.value(x))
    
    def __mul__(self, other):
        return Function(lambda x: self.value(x) * other.value(x))
    
    def __pow__(self, other):
        return Function(lambda x: self.value(x)**other)

    def derivative(self):
        return Function(lambda x: (self.f(x + self.dx) - self.f(x)) / self.dx)
    
    def integral(self, a, b):
        sum = 0
        x = a
        while x < b:
            sum += self.f(x) * self.dx
            x += self.dx
        return sum
    
    def tangent(self, x):
        k = self.derivative().value(x)
        m = self.value(x) - k * x
        return Function(lambda x: k * x + m)
    
    def zero(self):
        xs = []
        x = self.a + self.dx
        prev_y = self.value(self.a)
        while x <= self.b:
            y = self.value(x)
            if sign(prev_y) != sign(y) or y == 0:
                xs.append(x)
            x += self.dx
            prev_y = y
        return xs
    
    def intersect(self, other):
        return (self - other).zero()
    
    def extreme_points(self):
        return self.derivative().zero()


class Polynom:
    def __init__(self, pol) -> None:
        self.pol = pol # key = exponent, value = coefficient

    def to_string(self):
        result = ''
        s = dict(sorted(self.pol.items(), key=lambda item: item[0], reverse=True))
        for e, c in s.items():
            if c == 0: continue
            if len(result) > 0 or c < 0: result += ' + ' if c > 0 else ' - '
            if abs(c) != 1 or e == 0: result += str(abs(c))
            if e != 0: result += 'x' if e == 1 else f'x^{e}'
        return result
    
    def coef(self, exp): # if not sure if exp exists in dict
        if exp in self.pol: return self.pol[exp]
        return 0
    
    def is_null(self):
        return all(c == 0 for c in self.pol.values())

    def derivative(self):
        result = {}
        for e, c in self.pol.items():
            if e == 0: continue
            result[e - 1] = c * e
        return Polynom(result)

    def primitive(self):
        result = {}
        for e, c in self.pol.items():
            result[e + 1] = c / (e + 1)
        return Polynom(result)
    
    def value(self, x):
        return sum(self.pol[e] * x**e for e in self.pol)
    
    def integral(self, a, b):
        prim = self.primitive()
        return prim.value(b) - prim.value(a)
    
    def gcd(self, other):
        q, r = self / other
        if r.deg() == 0: return other
        return other.gcd(r)
    
    def zero(self):
        deg = self.deg()

        if deg == 0: return []
        if deg == 1:
            return [-self.pol[0] / self.pol[1]]
        if deg == 2:
            q, r = self / Polynom({0: self.pol[2]})
            a = -q.coef(1) / 2
            b = a**2 - q.coef(0)
            #if b < 0: return []
            b = b**0.5
            return [a + b, a - b]
        
        # heltalspolynom
        for c in self.pol.values():
            if c % 1 != 0:
                break
        else:
            ps = divisors(self.pol[0])
            qs = divisors(self.pol[deg])
            for p in ps:
                for q in qs:
                    x = p / q
                    if self.value(x) == 0: return [x] + (self / Polynom({1: 1, 0: -x}))[0].zero()
                    if self.value(-x) == 0: return [-x] + (self / Polynom({1: 1, 0: x}))[0].zero()
        
        
        for guess in range(-15, 15): # guesses
            if self.value(guess) == 0:
                return [guess] + (self / Polynom({1: 1, 0: -guess}))[0].zero()
        
        l = -15
        for r in range(-15, 15): # guesses
            ry = self.value(r)
            ly = self.value(l)
            if ly < 0 and ry > 0 or ry < 0 and ly > 0: break
            l = r
        else: return []
        while l < r:
            m = (l + r) / 2
            my = self.value(m)
            if abs(my) <= 0.000001: break
            ly = self.value(l)
            ry = self.value(r)
            if my > 0 and ry > 0 or my < 0 and ry < 0: r = m
            else: l = m
        return [m] + (self / Polynom({1: 1, 0: -m}))[0].zero()

        return [] # could not find

    def intersect(self, other):
        p = self - other
        return p.zero()

    def clean(self):
        result = {}
        for e, c in self.pol.items():
            if c != 0:
                result[e] = c
        self.pol = result
    
    def deg(self):
        record = 0
        for e, c in self.pol.items():
            if c != 0 and e > record: record = e
        return record

    def tangent(self, x):
        k = self.derivative().value(x)
        m = self.value(x) - k * x
        return Polynom({1: k, 0: m})
    
    def min(self):
        xs = self.derivative().zero()
        return min(xs, key=lambda x: self.value(x))
    
    def max(self):
        xs = self.derivative().zero()
        return max(xs, key=lambda x: self.value(x))
    
    def dist_to_pt(self, x, y):
        dist = Polynom({1: 1, 0: -x})**2 + (self - Polynom({0: y}))**2
        return dist.min() # x of min dist

    def __truediv__(self, other):
        q = Polynom({}) # kvot
        r = Polynom(self.pol.copy()) # rest
        other_deg = other.deg()
        r_deg = r.deg()
        while not r.is_null() and r_deg >= other_deg:
            deg = r_deg - other_deg
            coeff = r.coef(r_deg) / other.coef(other_deg)
            p = Polynom({deg: coeff})
            r -= p * other
            q += p
            r_deg = r.deg()
        return [q, r] # kvot, rest
    
    def __mul__(self, other):
        result = {}
        for e1, c1 in self.pol.items():
            for e2, c2 in other.pol.items():
                e = e1 + e2
                c = c1 * c2
                if e not in result: result[e] = 0
                result[e] += c
        return Polynom(result)
    
    def __add__(self, other):
        result = self.pol.copy()
        for e, c in other.pol.items():
            if e not in result: result[e] = 0
            result[e] += c
        return Polynom(result)
    
    def __sub__(self, other):
        result = self.pol.copy()
        for e, c in other.pol.items():
            if e not in result: result[e] = 0
            result[e] -= c
        return Polynom(result)
    
    def __pow__(self, other):
        result = Polynom({0: 1})
        for _ in range(other):
            result *= self
        return result
    


def sign(x):
    return -1 if x < 0 else 1

def line_fit(pts):
    mx = 0
    my = 0
    for pt in pts:
        mx += pt[0]
        my += pt[1]
    num_pts = len(pts)
    mx /= num_pts
    my /= num_pts
    
    dist_pol = Polynom({})
    for pt in pts:
        dist_pol += Polynom({1: mx - pt[0], 0: pt[1] - my})**2
    k = dist_pol.min()
    m = my - k * mx
    return Polynom({1: k, 0: m})


def divisors(n):
    return [i for i in range(1, int(abs(n) / 2) + 1) if n % i == 0] + [abs(n)]


def is_prime(x):
    for i in range(2, x // 2):
        if x % i == 0: return False
    return True


def factors(x):
    if is_prime(x):
        return [x]
    for i in range(2, x // 2):
        if x % i == 0:
            return factors(i) + factors(x // i)
