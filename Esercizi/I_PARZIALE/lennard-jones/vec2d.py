import numpy as np

## NOTE: definiamo la classe vettori in due dimensioni, la funzione _init_ serve a dare dei valori ad x ed
class vec2d:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    ## NOTE: possiamo fornire infiniti attributi ad una classe definendo delle funzioni
    def mod(self):
        return np.sqrt(self.x ** 2 + self.y ** 2)

    ##NOTE: funzione che definisce il versore corrispondente al vettore passato
    def unitary(self):
        modulo = self.mod()
        return vec2d(self.x / modulo, self.y / modulo)

    ##NOTE: funzione che definisce la somma tra vettori
    def __add__(self, other):
        return vec2d(self.x + other.x, self.y + other.y)

    ##NOTE: funzione che definisce la sottrazione tra vettori
    def __sub__(self, other):
        return vec2d(self.x - other.x, self.y - other.y)

    ##NOTE: funzione che definisce il prodotto tra uno scalare e un vettore o tra due vattori (scalare)
    def __mul__(self, other):

        if type(other) == vec2d:
            return self.x * other.x + self.y * other.y

        return vec2d(self.x * other, self.y * other)

    def __rmul__(self, other):
        return self * other

    ##NOTE: funzione che definisce la divisione di un vettore per uno scalare
    def __truediv__(self, other):
        return self * (1 / other)

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ')'

    ##NOTE: funzione che resituisce l'angolo tra due vettori
    def get_angle(self, other, unit):
        if unit == 'rad':
            return np.arccos(self * other / (self.mod() * other.mod()))

        elif unit == 'deg':
            return (np.arccos(self * other / (self.mod() * other.mod()))) * 180 / np.pi

        else:
            return False

    ##NOTE: funzione che restituisce la distanza tra due vettori 
    def get_distance(self, other):
        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)
        
    def abs_value(self):
        return vec2d(np.abs(self.x), np.abs(self.y))

    def __mod__(self, other):
        return vec2d(self.x % other, self.y % other)

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ')'
