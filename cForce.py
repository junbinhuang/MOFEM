import numpy as np 
import math

def customizedForce(x,y): # Boundary force
    return bracket(x,y)

def bracket(x,y):
    """Used in the AMORE bracket example."""
    assert(y<=1.0e-6),"Error in input!"
    p=2000.0
    x=x-27.0
    theta=abs(math.atan2(y,x)+math.pi/2)
    px=2*(math.pi/2-theta)/math.pi*p*x/math.sqrt(x**2+y**2)
    py=2*(math.pi/2-theta)/math.pi*p*y/math.sqrt(x**2+y**2)
    return np.array([[px],[py]])