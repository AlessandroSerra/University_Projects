import turtle

class MyTurtle(turtle.Turtle):
  def __init__(self, screen = turtle.Screen()):
    turtle.Turtle.__init__(self, screen)
    self.hideturtle()
    
def create_turtles(screen, n = 10):
  for i in range(n):
    MyTurtle(screen)
    
def move_turtles(screen, dist=10, angle = 4):
  for i, turtle in enumerate(screen.turtles()):
    turtle.left(angle*(1+i))
    turtle.forward(dist)
    x, y = turtle.pos()
    try:
      turtle.color(abs(x), abs(y), abs(x+y))
    except:
      pass
    
writer = MyTurtle()
writer.penup()
writer.goto(0,100)
writer.write("Click Me!", font=("Arial",30), align = "center")