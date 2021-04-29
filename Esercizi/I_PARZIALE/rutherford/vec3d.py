import numpy as np

## NOTE: definiamo la classe vettori in tre dimensioni, la funzione _init_ serve a dare dei valori ad x, y e z
class vec3d:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    ## NOTE: possiamo fornire infiniti attributi ad una classe definendo delle funzioni
    def mod(self):
        return np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def unitary(self):
        modulo = self.mod()
        return vec3d(self.x / modulo, self.y / modulo, self.z / modulo)

    def __add__(self, other):
        return vec3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return vec3d(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):

        if type(other) == vec3d:
            return self.x * other.x + self.y * other.y + self.z * other.z

        return vec3d(self.x * other, self.y * other, self.z * other)

    def __rmul__(self, other):
        return self * other

    def cross_prod(self, other):
        prod_x = self.y * other.z - self.z * other.y
        prod_y = self.z * other.x - self.x * other.z
        prod_z = self.x * other.y - self.y * other.x
        return vec3d(prod_x, prod_y, prod_z)

    def __truediv__(self, other):
        return self * (1 / other)

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'

    def get_angle(self, other, unit):
        if unit == 'rad':
            return np.arccos(self * other / (self.mod() * other.mod()))

        elif unit == 'deg':
            return (np.arccos(self * other / (self.mod() * other.mod()))) * 180 / np.pi

        else:
            print('Invalid unit, pick rad or deg\n')