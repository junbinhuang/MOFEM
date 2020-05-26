# Stiffness matrix calculation.
import numpy as np 
import scipy
import sys, os
sys.path.append(os.path.dirname(sys.path[0]))
from otherFunctions import numericalIntegration
from elementLibrary import shapeFunction,invShapeFunction

def triFE(coord,Dmat):
    nInt=1

    (pos,wei)=numericalIntegration.gaussTri(nInt)

    Kloc=np.zeros((2*len(coord),2*len(coord)))

    for i in range(nInt):
        (Bmat,Jacobian)=shapeFunction.shapeTriFE(coord,pos[i,0],pos[i,1])[1:3]
        Kloc+=0.5*np.matmul(Bmat.transpose(),np.matmul(Dmat,Bmat))*Jacobian*wei[i]
    
    return Kloc

def triQuadFE(coord,Dmat):
    nInt=3

    (pos,wei)=numericalIntegration.gaussTri(nInt)

    Kloc=np.zeros((2*len(coord),2*len(coord)))

    for i in range(nInt):
        (Bmat,Jacobian)=shapeFunction.shapeTriQuadFE(coord,pos[i,0],pos[i,1])[1:3]
        Kloc+=0.5*np.matmul(Bmat.transpose(),np.matmul(Dmat,Bmat))*Jacobian*wei[i]
    
    return Kloc

def quadFE(coord,Dmat):
    nInt=2 # 2 quadrature points in each direction.

    (pos,wei)=numericalIntegration.gaussQuad(nInt)

    Kloc=np.zeros((2*len(coord),2*len(coord)))

    for i in range(nInt):
        for j in range(nInt):
            _,Bmat,Jacobian=shapeFunction.shapeQuadFE(coord,pos[i],pos[j])
            Kloc+=np.matmul(Bmat.transpose(),np.matmul(Dmat,Bmat))*Jacobian*wei[i]*wei[j]
    
    return Kloc

def ICMFE(coord,Dmat): # For square elements
    nInt=2 # 2 quadrature points in each direction.

    (pos,wei)=numericalIntegration.gaussQuad(nInt)

    Kloc=np.zeros((2*len(coord),2*len(coord)))
    Hloc=np.zeros((4,4))
    Eloc=np.zeros((2*len(coord),4))

    for i in range(nInt):
        for j in range(nInt):
            Bmat,Gmat,Jacobian=shapeFunction.shapeICMFE(coord,pos[i],pos[j])
            Kloc+=np.matmul(Bmat.transpose(),np.matmul(Dmat,Bmat))*Jacobian*wei[i]*wei[j]
            Hloc+=np.matmul(Gmat.transpose(),np.matmul(Dmat,Gmat))*Jacobian*wei[i]*wei[j]
            Eloc+=np.matmul(Bmat.transpose(),np.matmul(Dmat,Gmat))*Jacobian*wei[i]*wei[j]

    Kloc-=Eloc@np.linalg.solve(Hloc,Eloc.transpose())

    return Kloc,Hloc,Eloc

def quadQuadFE(coord,Dmat):
    nInt=3 # 2 quadrature points in each direction.

    (pos,wei)=numericalIntegration.gaussQuad(nInt)

    Kloc=np.zeros((2*len(coord),2*len(coord)))

    for i in range(nInt):
        for j in range(nInt):
            _,Bmat,Jacobian=shapeFunction.shapeQuadQuadFE(coord,pos[i],pos[j])
            Kloc+=np.matmul(Bmat.transpose(),np.matmul(Dmat,Bmat))*Jacobian*wei[i]*wei[j]
    
    return Kloc

def sparsifyElementMatrix(K,numbering1,numbering2=None):
    """Convert the element stiffness matrix into COO format for assembling."""
    if numbering2 is None: numbering2=numbering1

    IArray=[]
    JArray=[]
    VArray=[]

    for i in range(2*len(numbering1)):
        for j in range(2*len(numbering2)):
            IArray.append(2*numbering1[i//2]+i%2)
            JArray.append(2*numbering2[j//2]+j%2)
            VArray.append(K[i,j])
    
    return IArray,JArray,VArray

def lowerOrderAMORE(coordTri,Dmat,coord1,rho1,coord2=None,rho2=None):
    if coord2 is None: 
        coord2=coord1
        rho2=rho1

    if coord1.shape[0]+coord2.shape[0]==6: nInt=3
    elif coord1.shape[0]+coord2.shape[0]==7: nInt=4
    else: nInt=6 

    pos,wei=numericalIntegration.gaussTri(nInt)
    LCmat1=getLC(coordTri,rho1)

    if rho2 is rho1: LCmat2=LCmat1
    else: LCmat2=getLC(coordTri,rho2)

    Kloc=np.zeros((2*coord1.shape[0],2*coord2.shape[0]))

    for i in range(nInt):
        Nmat,_,Jacobian=shapeFunction.shapeTriFE2(coordTri,pos[i,0],pos[i,1])
        rho1Value=(Nmat@rho1)[0]
        xy=(Nmat@coordTri).reshape(-1)
        Bmat1=getStrainMatrix(LCmat1,coord1,xy,rho1Value)

        if rho2 is rho1: Bmat2=Bmat1
        else: 
            rho2Value=(Nmat@rho2)[0]
            Bmat2=getStrainMatrix(LCmat2,coord2,xy,rho2Value)

        Kloc+=0.5*np.matmul(Bmat1.transpose(),np.matmul(Dmat,Bmat2))*Jacobian*wei[i]
    
    return Kloc

def quadraticAMORE(coordTri,Dmat,coord1,rho1,coord2=None,rho2=None):
    """For a straight-edge triangle. It can be curved only at the boundary."""
    if coord2 is None: 
        coord2=coord1
        rho2=rho1

    if coord1.shape[0]+coord2.shape[0]==6: nInt=3 # 3+3
    elif coord1.shape[0]+coord2.shape[0]==7: nInt=4 # 3+4
    elif coord1.shape[0]+coord2.shape[0]==8: nInt=6 # 4+4
    elif coord1.shape[0]+coord2.shape[0]==9: nInt=4 # 3+6
    elif coord1.shape[0]+coord2.shape[0]==10: nInt=6 # 4+6
    elif coord1.shape[0]+coord2.shape[0]==12: nInt=9 # 3+9
    elif coord1.shape[0]+coord2.shape[0]==13: nInt=12 # 4+9
    else: nInt=16 # 9+9

    pos,wei=numericalIntegration.gaussTri(nInt)

    Kloc=np.zeros((2*coord1.shape[0],2*coord2.shape[0]))

    if len(coordTri)==3:
        LCmat1=getLC(coordTri,rho1) # A constant matrix

        if rho2 is rho1: LCmat2=LCmat1
        else: LCmat2=getLC(coordTri,rho2)

        for i in range(nInt):
            Nmat,_,Jacobian=shapeFunction.shapeTriFE2(coordTri,pos[i,0],pos[i,1])
            rho1Value=(Nmat@rho1)[0]
            xy=(Nmat@coordTri).reshape(-1)
            Bmat1=getStrainMatrix(LCmat1,coord1,xy,rho1Value)
    
            if rho2 is rho1: Bmat2=Bmat1
            else: 
                rho2Value=(Nmat@rho2)[0]
                Bmat2=getStrainMatrix(LCmat2,coord2,xy,rho2Value)
    
            Kloc+=0.5*np.matmul(Bmat1.transpose(),np.matmul(Dmat,Bmat2))*Jacobian*wei[i]   

    elif len(coordTri)==6:
        for i in range(nInt):
            LCmat1=getLC(coordTri,rho1,pos[i,0],pos[i,1]) # No longer constant
            if rho2 is rho1: LCmat2=LCmat1
            else: LCmat2=getLC(coordTri,rho2,pos[i,0],pos[i,1])

            Nmat,_,Jacobian=shapeFunction.shapeTriQuadFE2(coordTri,pos[i,0],pos[i,1])
            rho1Value=pos[i,0]*rho1[0]+pos[i,1]*rho1[1]+pos[i,2]*rho1[2]

            xy=(Nmat@coordTri).reshape(-1)
            Bmat1=getStrainMatrix(LCmat1,coord1,xy,rho1Value)

            if rho2 is rho1: Bmat2=Bmat1
            else: 
                rho2Value=pos[i,0]*rho2[0]+pos[i,1]*rho2[1]+pos[i,2]*rho2[2]
                Bmat2=getStrainMatrix(LCmat2,coord2,xy,rho2Value)

            Kloc+=0.5*np.matmul(Bmat1.transpose(),np.matmul(Dmat,Bmat2))*Jacobian*wei[i]

    else: raise ValueError
    
    return Kloc

def getStrainMatrix(LCmat,coord,xy,rho):
    if coord.shape[0]==4: 
        isoCoord=invShapeFunction.invQuad(coord,xy)
        Nmat,Bmat,_=shapeFunction.shapeQuadFE(coord,isoCoord[0],isoCoord[1])
    elif coord.shape[0]==9:
        isoCoord=invShapeFunction.invQuadQuad(coord,xy)
        Nmat,Bmat,_=shapeFunction.shapeQuadQuadFE(coord,isoCoord[0],isoCoord[1])
    elif coord.shape[0]==6:
        isoCoord=invShapeFunction.invTriQuad(coord,xy)
        Nmat,Bmat,_=shapeFunction.shapeTriQuadFE(coord,isoCoord[0],isoCoord[1])     
    elif coord.shape[0]==3: 
        isoCoord=invShapeFunction.invTri(coord,xy)
        Nmat,Bmat,_=shapeFunction.shapeTriFE(coord,isoCoord[0],isoCoord[1])
    else: raise ValueError
    
    return (LCmat@Nmat) + (rho*Bmat)

def getLC(coordTri,rho,r=0.0,s=0.0):
    if len(coordTri)==3:
        dudx=shapeFunction.shapeTriFE2(coordTri,r,s)[1]
        drhodx=dudx@rho
    elif len(coordTri)==6:
        dudx=shapeFunction.shapeTriQuadFE2(coordTri,r,s)[1]
        rhoExtension=np.zeros(6)
        rhoExtension[0:3]=rho
        rhoExtension[3:6]=0.5*(rho+rho[[1,2,0]])
        drhodx=dudx@rhoExtension
    else: raise ValueError

    LC=np.array([[drhodx[0],0.0],[0.0,drhodx[1]],[drhodx[1],drhodx[0]]])
    return LC

def getCoordFromHalfEdge(edge):
    """Obtain the coordinates of a triangle in the half-edge form."""
    startingEdge=edge
    coord=[]
    coord.append([edge.origin.x,edge.origin.y])
    edge=edge.next

    while edge is not startingEdge:
        coord.append([edge.origin.x,edge.origin.y])
        edge=edge.next   
        
    coord=np.array(coord,dtype='d')
    return coord

if __name__=='__main__':
    pass