import numpy as np 
import math
import sys, os
sys.path.append(os.path.dirname(sys.path[0]))
from elementLibrary import shapeFunction

def invLine(coord,xy):
    """Input:
    
    coord=array([[x1,y1],[x2,y2],...]);
    
    xy=array([x,y]).
    
    xy must be between coord[0] and coord[1]."""

    length=lenEdge(coord[0],coord[1])
    length1=lenEdge(coord[0],xy)

    return 2.0*length1/length-1.0

def invTri(coord,xy):
    """Input:
    
    coord=array([[x1,y1],[x2,y2],[x3,y3]]);
    
    xy=array([x,y])."""

    area=(coord[1,0]*coord[2,1]-coord[2,0]*coord[1,1]-coord[0,0]*\
        (coord[2,1]-coord[1,1])+coord[0,1]*(coord[2,0]-coord[1,0]))/2

    a=np.array([coord[1,0]*coord[2,1]-coord[2,0]*coord[1,1],\
        coord[2,0]*coord[0,1]-coord[0,0]*coord[2,1],\
            coord[0,0]*coord[1,1]-coord[1,0]*coord[0,1]])

    b=np.array([coord[1,1]-coord[2,1],coord[2,1]-coord[0,1],coord[0,1]-coord[1,1]])

    c=np.array([coord[2,0]-coord[1,0],coord[0,0]-coord[2,0],coord[1,0]-coord[0,0]])

    return (a+b*xy[0]+c*xy[1])/2.0/area

def invQuad(coord,xy):
    """Input:
    
    coord=array([[x1,y1],[x2,y2],[x3,y3],[x4,y4]]);
    
    xy=array([x,y])."""

    elementSize=max(abs(coord[2,0]-coord[0,0]),abs(coord[2,1]-coord[0,1]))

    tol = elementSize * 10**(-15)
    isoCoord=np.zeros(2) # Initiation

    xyNew=xy-coord[0] # Normalization
    coordNew=coord.copy()
    for i in range(4):
        coordNew[i]-=coord[0]
    
    Nmat,dNmat=quadShapeFunction(isoCoord)

    count=0
    res=xyNew-(Nmat@coordNew).reshape(-1)
    while np.amax(np.absolute(res))>tol:
        J=(dNmat@coordNew).transpose()
        isoCoord+=np.linalg.solve(J,res)

        count+=1

        assert count<20, "Too many iterations!"

        Nmat,dNmat=quadShapeFunction(isoCoord)
        res=xyNew-(Nmat@coordNew).reshape(-1)
    
    return isoCoord

def invQuadQuad(coord,xy):
    """Input:
    
    coord=array([[x1,y1],[x2,y2],[x3,y3],[x4,y4],[x5,y5],[x6,y6],[x7,y7],[x8,y8],[x9,y9]]);
    
    xy=array([x,y])."""

    elementSize=max(abs(coord[2,0]-coord[0,0]),abs(coord[2,1]-coord[0,1]))

    tol = elementSize * 10**(-15)
    isoCoord=np.zeros(2) # Initiation

    xyNew=xy-coord[0] # Normalization
    coordNew=coord.copy()
    for i in range(len(coord)):
        coordNew[i]-=coord[0]
    
    Nmat,dNmat=quadQuadShapeFunction(isoCoord)

    count=0
    res=xyNew-(Nmat@coordNew).reshape(-1)
    while np.amax(np.absolute(res))>tol:
        J=(dNmat@coordNew).transpose()
        isoCoord+=np.linalg.solve(J,res)

        count+=1

        assert count<20, "Too many iterations!"

        Nmat,dNmat=quadQuadShapeFunction(isoCoord)
        res=xyNew-(Nmat@coordNew).reshape(-1)
    
    return isoCoord

def invTriQuad(coord,xy):
    """Input:
    
    coord=array([[x1,y1],[x2,y2],[x3,y3],[x4,y4],[x5,y5],[x6,y6]]);
    
    xy=array([x,y])."""

    elementSize=max(abs(coord[2,0]-coord[0,0]),abs(coord[2,1]-coord[0,1]))

    tol = elementSize * 10**(-15)
    isoCoord=np.array([1/3,1/3]) # Initiation

    xyNew=xy-coord[0] # Normalization
    coordNew=coord.copy()
    for i in range(len(coord)):
        coordNew[i]-=coord[0]
    
    Nmat,dNmat=triQuadShapeFunction(isoCoord)

    count=0
    res=xyNew-(Nmat@coordNew).reshape(-1)
    while np.amax(np.absolute(res))>tol:
        J=(dNmat@coordNew).transpose()
        isoCoord+=np.linalg.solve(J,res)

        count+=1

        assert count<20, "Too many iterations!"

        Nmat,dNmat=triQuadShapeFunction(isoCoord)
        res=xyNew-(Nmat@coordNew).reshape(-1)
    
    return isoCoord

def quadShapeFunction(isoCoord):
    """Input: 
    
    isoCoord=array([r,s])."""

    Nmat=np.array([[0.25*(1.0-isoCoord[0])*(1.0-isoCoord[1]),\
                    0.25*(1.0+isoCoord[0])*(1.0-isoCoord[1]),\
                    0.25*(1.0+isoCoord[0])*(1.0+isoCoord[1]),\
                    0.25*(1.0-isoCoord[0])*(1.0+isoCoord[1])]])

    dNmat=np.array([[-0.25*(1.0-isoCoord[1]),0.25*(1.0-isoCoord[1]),\
                    0.25*(1.0+isoCoord[1]),-0.25*(1.0+isoCoord[1])],
                    [-0.25*(1.0-isoCoord[0]),-0.25*(1.0+isoCoord[0]),\
                    0.25*(1.0+isoCoord[0]),0.25*(1.0-isoCoord[0])]])

    return Nmat,dNmat

def quadQuadShapeFunction(isoCoord):
    """Input: 
    
    isoCoord=array([r,s])."""
    r=isoCoord[0]
    s=isoCoord[1]

    Nmat=np.array([[0.25*(1.0-r)*(1.0-s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0-s)*(1.0-r**2)-0.25*(1.0-r)*(1.0-s**2),\
                    0.25*(1.0+r)*(1.0-s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0-s)*(1.0-r**2)-0.25*(1.0+r)*(1.0-s**2),\
                    0.25*(1.0+r)*(1.0+s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0+r)*(1.0-s**2)-0.25*(1.0+s)*(1.0-r**2),\
                    0.25*(1.0-r)*(1.0+s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0+s)*(1.0-r**2)-0.25*(1.0-r)*(1.0-s**2),\
                    0.5*(1.0-s)*(1.0-r**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    0.5*(1.0+r)*(1.0-s**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    0.5*(1.0+s)*(1.0-r**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    0.5*(1.0-r)*(1.0-s**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    (1.0-r**2)*(1.0-s**2)]])

    dNmat=np.array([[-0.25*(1.0-s)-0.5*r*(1.0-s**2)+0.5*r*(1.0-s)+0.25*(1.0-s**2),\
                     0.25*(1.0-s)-0.5*r*(1.0-s**2)+0.5*r*(1.0-s)-0.25*(1.0-s**2),\
                     0.25*(1.0+s)-0.5*r*(1.0-s**2)-0.25*(1.0-s**2)+0.5*r*(1.0+s),\
                     -0.25*(1.0+s)-0.5*r*(1.0-s**2)+0.5*r*(1.0+s)+0.25*(1.0-s**2),\
                     -r*(1.0-s)+r*(1.0-s**2),\
                     0.5*(1.0-s**2)+r*(1.0-s**2),\
                     -r*(1.0+s)+r*(1.0-s**2),\
                     -0.5*(1.0-s**2)+r*(1.0-s**2),\
                     -2.0*r*(1.0-s**2)],\
                    [-0.25*(1.0-r)-0.5*s*(1.0-r**2)+0.25*(1.0-r**2)+0.5*s*(1.0-r),\
                     -0.25*(1.0+r)-0.5*s*(1.0-r**2)+0.25*(1.0-r**2)+0.5*s*(1.0+r),\
                     0.25*(1.0+r)-0.5*s*(1.0-r**2)+0.5*s*(1.0+r)-0.25*(1.0-r**2),\
                     0.25*(1.0-r)-0.5*s*(1.0-r**2)-0.25*(1.0-r**2)+0.5*s*(1.0-r),\
                     -0.5*(1.0-r**2)+s*(1.0-r**2),\
                     -s*(1.0+r)+s*(1.0-r**2),\
                     0.5*(1.0-r**2)+s*(1.0-r**2),\
                     -s*(1.0-r)+s*(1.0-r**2),\
                     -2.0*s*(1.0-r**2)]])

    return Nmat,dNmat

def triShapeFunction(isoCoord):
    """Input: 
    
    isoCoord=array([L1,L2])."""

    Nmat=np.array([[isoCoord[0],isoCoord[1],1.0-isoCoord[0]-isoCoord[1]]])

    dNmat=np.array([[1.0,0.0,-1.0],
                    [0.0,1.0,-1.0]])

    return Nmat,dNmat

def triQuadShapeFunction(isoCoord):
    """Input: 
    
    isoCoord=array([L1,L2])."""
    r=isoCoord[0]
    s=isoCoord[1]
    t=1.0-r-s

    Nmat=np.array([[r-2.0*r*s-2.0*r*t,\
                    s-2.0*s*t-2.0*r*s,\
                    t-2.0*s*t-2.0*r*t,\
                    4.0*r*s,4.0*s*t,4.0*r*t]])

    dNmat=np.array([[1.0-2.0*s-2.0*t+2.0*r,\
                     0.0,\
                     -1.0+2.0*s-2.0*t+2.0*r,\
                     4.0*s,-4.0*s,4.0*t-4.0*r],\
                    [0.0,\
                     1.0-2.0*t+2.0*s-2.0*r,\
                     -1.0-2.0*t+2.0*s+2.0*r,\
                     4.0*r,4.0*t-4.0*s,-4.0*r]])

    return Nmat,dNmat

def lenEdge(coord1,coord2):
    """Input:
    
    coord1=array([x1,y1]);
    
    coord2=array([x2,y2])."""

    return math.sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2)

if __name__=='__main__':
    pass