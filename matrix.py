


class Mat:
    def __init__(self, m):
        self.m = m

    def size(self):
        return len(self.m), len(self.m[0]) # rows, cols
    
    def transpose(self):
        rs, cs = self.size()
        return Mat([[self.m[r][c]  for r in range(rs)] for c in range(cs)])
    
    def col(self, i):
        rs, _ = self.size()
        return [self.m[r][i] for r in range(rs)]
    
    def copy(self):
        return Mat([[x for x in r] for r in self.m])
    
    def flip(self):
        return Mat([r[::-1] for r in self.m[::-1]])
    
    @staticmethod
    def unit(n):
        return Mat([[int(i == j) for j in range(n)] for i in range(n)])

    @staticmethod
    def all(v, n, m):
        return Mat([[v for _ in range(m)] for _ in range(n)])

    @staticmethod
    def zeros(n, m):
        return Mat.all(0, n, m)

    @staticmethod
    def ones(n, m):
        return Mat.all(1, n, m)
    
    def __repr__(self):
        s = ''
        for r in self.m:
            s += f"[{' '.join(str(x) for x in r)}]\n"
        return s[:-1]
    
    def __add__(self, other):
        ar, ac = self.size()
        br, bc = other.size()
        if ar != br or ac != bc:
            raise ValueError('size(a) != size(b)')
        return Mat([[self.m[r][c] + other.m[r][c] for c in range(ac)] for r in range(ar)])

    def __sub__(self, other):
        ar, ac = self.size()
        br, bc = other.size()
        if ar != br or ac != bc:
            raise ValueError('size(a) != size(b)')
        return Mat([[self.m[r][c] - other.m[r][c] for c in range(ac)] for r in range(ar)])
    
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Mat([[x * other for x in r] for r in self.m])
        ar, ac = self.size()
        br, bc = other.size()
        if ac != br:
            raise ValueError('cols(a) != rows(b)')
        return Mat([[dot(self.m[r], other.col(c)) for c in range(bc)] for r in range(ar)])

    def __rmul__(self, other):
        return self * other

    def __pow__(self, other):
        if not isinstance(other, int):
            raise ValueError('exponent must be int')
        if other < 0:
            raise ValueError('negative exponent')
        rs, cs = self.size()
        if rs != cs:
            raise ValueError('Rows != Cols')
        if other == 0:
            return Mat.unit(rs)
        m = self.copy()
        for _ in range(other - 1):
            m *= self
        return m

    def invert(self): # gauss
        m = self.copy()
        rs, cs = m.size()
        if rs != cs:
            raise ValueError('Rows != Cols')
        n = rs
        b = Mat.unit(n)
        m, b = gauss(m, b)
        m, b = gauss(m.flip(), b.flip())
        return b.flip()

class ColVec(Mat):
    def __init__(self, *v: float):
        super().__init__([[x] for x in v])

    def __repr__(self):
        s = ''
        for r in self.m:
            s += f"[{r[0]}]\n"
        return s[:-1]

class RowVec(Mat):
    def __init__(self, *v: float):
        super().__init__([v])

    def __repr__(self):
        s = ''
        for r in self.m:
            s += f"[{' '.join(str(x) for x in r)}]\n"
        return s[:-1]


def dot(a, b):
    return sum(i * j for i, j in zip(a, b))



def gauss(m, b):
    m = m.copy()
    b = b.copy()
    mr, mc = m.size()
    br, bc = b.size()
    if mr != br:
        raise ValueError('mr != br')
    n = min(mr, mc)
    for i in range(n):
        for j in range(i + 1, n):
            s = -m.m[j][i] / m.m[i][i] # scalar
            for k in range(mc):
                m.m[j][k] += m.m[i][k] * s
            for k in range(bc):
                b.m[j][k] += b.m[i][k] * s
    for i in range(n):
        s = m.m[i][i]
        for k in range(mc):
            m.m[i][k] /= s
        for k in range(bc):
            b.m[i][k] /= s
    return m, b


