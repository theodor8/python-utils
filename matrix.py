


class Matrix:
    def __init__(self, m):
        self.m = m

    def size(self):
        return len(self.m), len(self.m[0]) # rows, cols
    
    def transpose(self):
        rs, cs = self.size()
        return Matrix([[self.m[r][c]  for r in range(rs)] for c in range(cs)])
    
    def col(self, i):
        rs, cs = self.size()
        return [self.m[r][i] for r in range(rs)]
    
    def copy(self):
        return Matrix([[x for x in r] for r in self.m])
    
    def flip(self):
        return Matrix([r[::-1] for r in self.m[::-1]])
    
    @staticmethod
    def unit(n):
        return Matrix([[int(i == j)  for j in range(n)] for i in range(n)])
    
    def __str__(self):
        s = ''
        for r in self.m:
            s += f"[{' '.join(str(x) for x in r)}]\n"
        return s[:-1]
    
    def __add__(self, other):
        ar, ac = self.size()
        br, bc = other.size()
        if ar != br or ac != bc:
            raise ValueError('size(a) != size(b)')
        return [[self.m[r][c] + other.m[r][c] for c in range(ac)] for r in range(ar)]
    
    def __mul__(self, other):
        ar, ac = self.size()
        br, bc = other.size()
        if ac != br:
            raise ValueError('cols(a) != rows(b)')
        return Matrix([[dot(self.m[r], other.col(c)) for c in range(bc)] for r in range(ar)])
    
    def invert(self): # gauss
        m = self.copy()
        rs, cs = m.size()
        if rs != cs:
            raise ValueError('Rows != Cols')
        n = rs
        b = Matrix.unit(n)
        m, b = gauss(m, b)
        m, b = gauss(m.flip(), b.flip())
        return b.flip()



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



m = Matrix([[-2, -3, 2, -3],
            [1, 0, -3, -2],
            [0, 1, -1, -1],
            [-2, 2, -1, -3]])

mi = m.invert()
b = Matrix([[1],
            [2],
            [3],
            [4]])
x = mi * b
print(m * x)
