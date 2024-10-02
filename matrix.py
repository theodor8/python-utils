

def size(m):
    return len(m), len(m[0]) # rows, cols

def dot(a, b):
    return sum(i * j for i, j in zip(a, b))

def transpose(m):
    rs, cs = size(m)
    return [[m[r][c]  for r in range(rs)] for c in range(cs)]

def add(a, b):
    ar, ac = size(a)
    br, bc = size(b)
    if ar != br or ac != bc:
        raise ValueError('size(a) != size(b)')
    return [[a[r][c] + b[r][c] for c in range(ac)] for r in range(ar)]

def scalar_mult(m, s):
    return [[x * s for x in r] for r in m]

def col(m, i):
    rs, cs = size(m)
    return [m[r][i] for r in range(rs)]

def unit_matrix(n):
    return [[int(i == j)  for j in range(n)] for i in range(n)]

def print_matrix(m):
    for r in m:
        print(f"[{' '.join(str(x) for x in r)}]")

def matmul(a, b):
    ar, ac = size(a)
    br, bc = size(b)
    if ac != br:
        raise ValueError('cols(a) != rows(b)')
    return [[dot(a[r], col(b, c)) for c in range(bc)] for r in range(ar)]

def matrix_copy(m):
    return [[x for x in r] for r in m]

def flip(m):
    return [r[::-1] for r in m[::-1]]

def gauss(m, b):
    m = matrix_copy(m)
    b = matrix_copy(b)
    mr, mc = size(m)
    br, bc = size(b)
    if mr != br:
        raise ValueError('mr != br')
    n = min(mr, mc)

    for i in range(n):
        for j in range(i + 1, n):
            s = -m[j][i] / m[i][i] # scalar
            for k in range(mc):
                m[j][k] += m[i][k] * s
            for k in range(bc):
                b[j][k] += b[i][k] * s
    for i in range(n):
        s = m[i][i]
        for k in range(mc):
            m[i][k] /= s
        for k in range(bc):
            b[i][k] /= s
    return m, b


def invert_matrix(m): # gauss
    m = matrix_copy(m)
    rs, cs = size(m)
    if rs != cs:
        raise ValueError('Rows != Cols')
    n = rs
    b = unit_matrix(n)
    m, b = gauss(m, b)
    m, b = gauss(flip(m), flip(b))
    return flip(b)



m = [[1, 2, 3, -1],
     [5, 6, 7, -2]]
b = [[6], [11]]

m, b = gauss(m, b)

print_matrix(m)
print_matrix(b)

for w in range(-100, 100):
    for z in range(-100, 100):
        y = 4.75 + 0.75*w - 2*z
        x = 6 + w - 3*z - 2*y
        if x + 2*y + 3*z - w != 6 or  5*x + 6*y + 7*z - 2*w != 11:
            print(1)