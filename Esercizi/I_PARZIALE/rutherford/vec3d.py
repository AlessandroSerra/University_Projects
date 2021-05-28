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

    ##NOTE: funzione che definisce il versore corrispondente al vettore passato
    def unitary(self):
        modulo = self.mod()
        return vec3d(self.x / modulo, self.y / modulo, self.z / modulo)

    ##NOTE: funzione che definisce l'opposto di un vettore
    def __neg__(self):
        return vec3d(-self.x, -self.y, -self.z)

    ##NOTE: funzione che definisce la somma tra vettori
    def __add__(self, other):
        return vec3d(self.x + other.x, self.y + other.y, self.z + other.z)

    ##NOTE: funzione che definisce la sottrazione tra vettori
    def __sub__(self, other):
        return vec3d(self.x - other.x, self.y - other.y, self.z - other.z)

    ##NOTE: funzione che definisce il prodotto tra uno scalare e un vettore o tra due vattori (scalare)
    def __mul__(self, other):

        if type(other) == vec3d:
            return self.x * other.x + self.y * other.y + self.z * other.z

        return vec3d(self.x * other, self.y * other, self.z * other)

    def __rmul__(self, other):
        return self * other

    ##NOTE: funzione che esegue il prodotto vettoriale mediante la definizione del tensore di Levi-Civita
    def cross_prod(self, other):
        prod_x = self.y * other.z - self.z * other.y
        prod_y = self.z * other.x - self.x * other.z
        prod_z = self.x * other.y - self.y * other.x
        return vec3d(prod_x, prod_y, prod_z)

    ##NOTE: funzione che definisce la divisione di un vettore per uno scalare
    def __truediv__(self, other):
        return self * (1 / other)

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'

    ##NOTE: funzione che restituisce la distanza tra due vettori 
    def get_distance(self, other):
        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other)**2)
    
    ##NOTE: funzione che restituisce il valore assoluto delle componenti di un vettore
    def abs_value(self):
        return vec3d(np.abs(self.x), np.abs(self.y), np.abs(self.z))

    ##NOTE: funzione che resituisce l'angolo tra due vettori
    def get_angle(self, other, unit):
        if unit == 'rad':
            return np.arccos(self * other / (self.mod() * other.mod()))

        elif unit == 'deg':
            return (np.arccos(self * other / (self.mod() * other.mod()))) * 180 / np.pi

        else:
            return False