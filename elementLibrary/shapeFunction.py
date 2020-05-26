# Shape functions for finite elements.
import numpy as np 
import math

def oneDLinear(r):
    Nmat=np.array([[0.5*(1.0-r),0.0,0.5*(1.0+r),0.0],\
        [0.0,0.5*(1.0-r),0.0,0.5*(1.0+r)]])

    return Nmat

def oneDLinear2(r):
    Nmat=np.array([[0.5*(1.0-r),0.5*(1.0+r)]])

    return Nmat

def oneDQuadratic(r,coord):
    Nmat=np.array([[0.5*(1.0-r)-0.5*(1.0-r**2),0.0,0.5*(1.0+r)-0.5*(1.0-r**2),0.0,1.0-r**2,0.0],\
                    [0.0,0.5*(1.0-r)-0.5*(1.0-r**2),0.0,0.5*(1.0+r)-0.5*(1.0-r**2),0.0,1.0-r**2]])
    dNmat=np.array([[-0.5+r,0.5+r,-2.0*r]])
    vector=dNmat@coord

    Jacobian=(vector[0,0]**2+vector[0,1]**2)**0.5

    return Nmat,Jacobian

def oneDQuadratic2(r):
    Nmat=np.array([[0.5*(1.0-r)-0.5*(1.0-r**2),0.5*(1.0+r)-0.5*(1.0-r**2),1.0-r**2]])

    return Nmat

def shapeQuadFE(coord,r,s):
    """Shape function of 4-node quadrilateral element."""

    Nmat=np.array([[0.25*(1.0-r)*(1.0-s),0.0,0.25*(1.0+r)*(1.0-s),0.0,\
        0.25*(1.0+r)*(1.0+s),0.0,0.25*(1.0-r)*(1.0+s),0.0],\
            [0.0,0.25*(1.0-r)*(1.0-s),0.0,0.25*(1.0+r)*(1.0-s),\
        0.0,0.25*(1.0+r)*(1.0+s),0.0,0.25*(1.0-r)*(1.0+s)]])

    dNmat=np.array([[-0.25*(1.0-s),0.25*(1.0-s),0.25*(1.0+s),-0.25*(1.0+s)],\
        [-0.25*(1.0-r),-0.25*(1.0+r),0.25*(1.0+r),0.25*(1.0-r)]])

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)

    Bmat=np.array([[dudx[0,0],0.0,dudx[0,1],0.0,dudx[0,2],0.0,dudx[0,3],0.0],\
        [0.0,dudx[1,0],0.0,dudx[1,1],0.0,dudx[1,2],0.0,dudx[1,3]],\
            [dudx[1,0],dudx[0,0],dudx[1,1],dudx[0,1],dudx[1,2],dudx[0,2],dudx[1,3],dudx[0,3]]])
    
    return Nmat,Bmat,Jacobian

def shapeICMFE(coord,r,s):
    """Shape function of 4-node quadrilateral element."""
    dNmat=np.array([[-0.25*(1.0-s),0.25*(1.0-s),0.25*(1.0+s),-0.25*(1.0+s)],\
        [-0.25*(1.0-r),-0.25*(1.0+r),0.25*(1.0+r),0.25*(1.0-r)]])

    dGmat=np.array([[-2.0*r,0.0],[0.0,-2.0*s]]) # No modification because we assume the element is a square.

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)
    dGdx=np.linalg.solve(J,dGmat)

    Bmat=np.array([[dudx[0,0],0.0,dudx[0,1],0.0,dudx[0,2],0.0,dudx[0,3],0.0],\
        [0.0,dudx[1,0],0.0,dudx[1,1],0.0,dudx[1,2],0.0,dudx[1,3]],\
            [dudx[1,0],dudx[0,0],dudx[1,1],dudx[0,1],dudx[1,2],dudx[0,2],dudx[1,3],dudx[0,3]]])

    Gmat=np.array([[dGdx[0,0],0.0,dGdx[0,1],0.0],\
                    [0.0,dGdx[1,0],0.0,dGdx[1,1]],\
                    [dGdx[1,0],dGdx[0,0],dGdx[1,1],dGdx[0,1]]])
    
    return Bmat,Gmat,Jacobian

def shapeQuadFE2(coord,r,s):
    """Shape function of 4-node quadrilateral element."""

    Nmat=np.array([[0.25*(1.0-r)*(1.0-s),0.25*(1.0+r)*(1.0-s),\
        0.25*(1.0+r)*(1.0+s),0.25*(1.0-r)*(1.0+s)]])

    dNmat=np.array([[-0.25*(1.0-s),0.25*(1.0-s),0.25*(1.0+s),-0.25*(1.0+s)],\
        [-0.25*(1.0-r),-0.25*(1.0+r),0.25*(1.0+r),0.25*(1.0-r)]])

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)
    
    return Nmat,dudx,Jacobian

def shapeQuadQuadFE(coord,r,s):
    """Shape function of 9-node quadrilateral element."""

    Nmat=np.array([[0.25*(1.0-r)*(1.0-s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0-s)*(1.0-r**2)-0.25*(1.0-r)*(1.0-s**2),0.0,\
                    0.25*(1.0+r)*(1.0-s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0-s)*(1.0-r**2)-0.25*(1.0+r)*(1.0-s**2),0.0,\
                    0.25*(1.0+r)*(1.0+s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0+r)*(1.0-s**2)-0.25*(1.0+s)*(1.0-r**2),0.0,\
                    0.25*(1.0-r)*(1.0+s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0+s)*(1.0-r**2)-0.25*(1.0-r)*(1.0-s**2),0.0,\
                    0.5*(1.0-s)*(1.0-r**2)-0.5*(1.0-r**2)*(1.0-s**2),0.0,\
                    0.5*(1.0+r)*(1.0-s**2)-0.5*(1.0-r**2)*(1.0-s**2),0.0,\
                    0.5*(1.0+s)*(1.0-r**2)-0.5*(1.0-r**2)*(1.0-s**2),0.0,\
                    0.5*(1.0-r)*(1.0-s**2)-0.5*(1.0-r**2)*(1.0-s**2),0.0,\
                    (1.0-r**2)*(1.0-s**2),0.0],\
                   [0.0,0.25*(1.0-r)*(1.0-s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0-s)*(1.0-r**2)-0.25*(1.0-r)*(1.0-s**2),\
                    0.0,0.25*(1.0+r)*(1.0-s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0-s)*(1.0-r**2)-0.25*(1.0+r)*(1.0-s**2),\
                    0.0,0.25*(1.0+r)*(1.0+s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0+r)*(1.0-s**2)-0.25*(1.0+s)*(1.0-r**2),\
                    0.0,0.25*(1.0-r)*(1.0+s)+0.25*(1.0-r**2)*(1.0-s**2)-0.25*(1.0+s)*(1.0-r**2)-0.25*(1.0-r)*(1.0-s**2),\
                    0.0,0.5*(1.0-s)*(1.0-r**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    0.0,0.5*(1.0+r)*(1.0-s**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    0.0,0.5*(1.0+s)*(1.0-r**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    0.0,0.5*(1.0-r)*(1.0-s**2)-0.5*(1.0-r**2)*(1.0-s**2),\
                    0.0,(1.0-r**2)*(1.0-s**2)]])

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

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)

    Bmat=np.array([[dudx[0,0],0.0,dudx[0,1],0.0,dudx[0,2],0.0,dudx[0,3],0.0,\
                    dudx[0,4],0.0,dudx[0,5],0.0,dudx[0,6],0.0,dudx[0,7],0.0,dudx[0,8],0.0],\
                   [0.0,dudx[1,0],0.0,dudx[1,1],0.0,dudx[1,2],0.0,dudx[1,3],\
                    0.0,dudx[1,4],0.0,dudx[1,5],0.0,dudx[1,6],0.0,dudx[1,7],0.0,dudx[1,8]],\
                   [dudx[1,0],dudx[0,0],dudx[1,1],dudx[0,1],dudx[1,2],dudx[0,2],dudx[1,3],dudx[0,3],\
                    dudx[1,4],dudx[0,4],dudx[1,5],dudx[0,5],dudx[1,6],dudx[0,6],dudx[1,7],dudx[0,7],dudx[1,8],dudx[0,8]]])
    
    return Nmat,Bmat,Jacobian

def shapeQuadQuadFE2(coord,r,s):
    """Shape function of 9-node quadrilateral element."""

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

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)
    
    return Nmat,dudx,Jacobian

def shapeTriFE(coord,r,s):
    """Shape function of 3-node triangular element."""

    Nmat=np.array([[r,0.0,s,0.0,1.0-r-s,0.0],\
        [0.0,r,0.0,s,0.0,1.0-r-s]])

    dNmat=np.array([[1.0,0.0,-1.0],\
        [0.0,1.0,-1.0]])

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)

    Bmat=np.array([[dudx[0,0],0.0,dudx[0,1],0.0,dudx[0,2],0.0],\
                   [0.0,dudx[1,0],0.0,dudx[1,1],0.0,dudx[1,2]],\
                   [dudx[1,0],dudx[0,0],dudx[1,1],dudx[0,1],dudx[1,2],dudx[0,2]]])

    return Nmat,Bmat,Jacobian

def shapeTriFE2(coord,r,s):
    """Shape function of 3-node triangular element."""

    Nmat=np.array([[r,s,1.0-r-s]])

    dNmat=np.array([[1.0,0.0,-1.0],\
        [0.0,1.0,-1.0]])

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)

    return Nmat,dudx,Jacobian

def shapeTriQuadFE(coord,r,s):
    """Shape function of 6-node triangular element."""
    t=1.0-r-s

    Nmat=np.array([[r-2.0*r*s-2.0*r*t,0.0,\
                    s-2.0*s*t-2.0*r*s,0.0,\
                    t-2.0*s*t-2.0*r*t,0.0,\
                    4.0*r*s,0.0,4.0*s*t,0.0,4.0*r*t,0.0],\
                   [0.0,r-2.0*r*s-2.0*r*t,\
                    0.0,s-2.0*s*t-2.0*r*s,\
                    0.0,t-2.0*s*t-2.0*r*t,\
                    0.0,4.0*r*s,0.0,4.0*s*t,0.0,4.0*r*t]])

    dNmat=np.array([[1.0-2.0*s-2.0*t+2.0*r,\
                     0.0,\
                     -1.0+2.0*s-2.0*t+2.0*r,\
                     4.0*s,-4.0*s,4.0*t-4.0*r],\
                    [0.0,\
                     1.0-2.0*t+2.0*s-2.0*r,\
                     -1.0-2.0*t+2.0*s+2.0*r,\
                     4.0*r,4.0*t-4.0*s,-4.0*r]])

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)

    Bmat=np.array([[dudx[0,0],0.0,dudx[0,1],0.0,dudx[0,2],0.0,dudx[0,3],0.0,dudx[0,4],0.0,dudx[0,5],0.0],\
                   [0.0,dudx[1,0],0.0,dudx[1,1],0.0,dudx[1,2],0.0,dudx[1,3],0.0,dudx[1,4],0.0,dudx[1,5]],\
                   [dudx[1,0],dudx[0,0],dudx[1,1],dudx[0,1],dudx[1,2],dudx[0,2],\
                    dudx[1,3],dudx[0,3],dudx[1,4],dudx[0,4],dudx[1,5],dudx[0,5]]])

    return Nmat,Bmat,Jacobian

def shapeTriQuadFE2(coord,r,s):
    """Shape function of 6-node triangular element."""
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

    J=np.matmul(dNmat,coord)
    Jacobian=np.linalg.det(J)

    dudx=np.linalg.solve(J,dNmat)

    return Nmat,dudx,Jacobian

if __name__=='__main__':
    pass