class Property:
    def __init__(self, x, y, type=None):
        self.x = x
        self.y = y
        self.type = type

    def test(self):
        return (self.x, self.y)

class City(Property):
    def test(self):
        return 'Hello'


p1 = Property(1, 2, 'blue')
print(p1.test())
c1 = City(p1)
print(c1.test())