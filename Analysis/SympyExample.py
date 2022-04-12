import sympy as sp
from PIL._imaging import display

#sp.init_printing(use_latex="mathjax")
sp.init_printing(2+sp.sqrt(3))

x=sp.symbols("x")
eq=sp.Eq(x**2-4,0)
display("x",eq)