


class Vector:
    def __init__(self, m): # m = array of values
        super().__init__([[v] for v in m])



    # def __rmul__(self, other):
    #     return self * other

    # def __truediv__(self, other):
    #     return Vector(*(v / other for v in self.vs))

    # def __add__(self, other):
    #     return Vector(*(i + j for i, j in zip(self.vs, other.vs)))

    # def __sub__(self, other):
    #     return Vector(*(i - j for i, j in zip(self.vs, other.vs)))

    # def __neg__(self):
    #     return Vector(*(-v for v in self.vs))
    
    # def cross(self, other):  # for 3 dims
    #     return Vector(self.vs[1] * other.vs[2] - self.vs[2] * other.vs[1],
    #                   self.vs[2] * other.vs[0] - self.vs[0] * other.vs[2],
    #                   self.vs[0] * other.vs[1] - self.vs[1] * other.vs[0])

    # def len(self): 
    #     return sum(r[0]**2 for r in self.m)**0.5

    # def proj(self, other):
    #     return ((other * self) / other.len()**2) * other

    # def angle(self, other):
    #     return math.acos((self * other) / (self.len() * other.len()))

    # def norm(self):
    #     return self / self.len()
    



