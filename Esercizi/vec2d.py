import numpy as np

## NOTE: definiamo la classe vettori in due dimensioni, la funzione _init_ serve a dare dei valori ad x ed
class vec2d:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    ## NOTE: possiamo fornire infiniti attributi ad una classe definendo delle funziono
    def mod(self):
        return np.sqrt(self.x ** 2 + self.y ** 2)

    def unitary(self):
        modulo = self.mod()
        return vec2d(self.x / modulo, self.y / modulo)

    def __add__(self, other):
        return vec2d(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return vec2d(self.x - other.x, self.y - other.y)

    def __mul__(self, other):

        if type(other) == vec2d:
            return self.x * other.x + self.y * other.y

        return vec2d(self.x * other, self.y * other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        return self * (1 / other)

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ')'

    def get_angle(self, other):
        return np.arccos(self * other / (self.mod() * other.mod()))