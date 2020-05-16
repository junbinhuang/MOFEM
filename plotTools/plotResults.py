import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import math
from matplotlib import rc
import sys, os
sys.path.append(os.path.dirname(sys.path[0]))
from linearSolvers import AMORE
from elementLibrary import shapeFunction, invShapeFunction, stiffnessMatrix
from plotTools import pM
from linearSolvers import AMORE

def getGraphDataOFE(coord,coordinates,incidentElements,edge,displacement,DMatrix,nSampling):
    X=np.zeros((nSampling,nSampling))
    Y=np.zeros((nSampling,nSampling))
    U=np.zeros((nSampling,nSampling))
    V=np.zeros((nSampling,nSampling))
    Sxx=np.zeros((nSampling,nSampling))
    Syy=np.zeros((nSampling,nSampling))
    Sxy=np.zeros((nSampling,nSampling))

    r=np.linspace(-1,1,nSampling)
    s=np.linspace(-1,1,nSampling)

    if len(coord)==3:
        newCoord=np.zeros((4,2))
        newCoord[0:3,:]=coord
        newCoord[3,:]=coord[2,:]

        for i in range(nSampling):
            for j in range(nSampling):
                Nmat,_=invShapeFunction.quadShapeFunction([r[i],s[j]])
                xy=(Nmat@newCoord).reshape(-1)
                X[i,j]=xy[0]
                Y[i,j]=xy[1]

                isoCoord=invShapeFunction.invTri(coord,xy)

                for k in range(len(incidentElements)):
                    if incidentElements[k]:
                        numbering=incidentElements[k].numbering
                        fNumbering=getFullNumber(numbering)
                        tempCoord=AMORE.getCoord(coordinates,numbering)
                        rho=AMORE.getRho(edge,k)
                        LCmat=stiffnessMatrix.getLC(coord,rho)
                        rhoValue=isoCoord[0]*rho[0]+isoCoord[1]*rho[1]+(1.0-isoCoord[0]-isoCoord[1])*rho[2]

                        if len(tempCoord)==4:
                            newIsoCoord=invShapeFunction.invQuad(tempCoord,xy)
                            Nmat,Bmat,_=shapeFunction.shapeQuadFE(tempCoord,newIsoCoord[0],newIsoCoord[1])
                            localDisp=rhoValue*Nmat@displacement[fNumbering]
                            localStress=rhoValue*DMatrix@Bmat@displacement[fNumbering]+DMatrix@LCmat@Nmat@displacement[fNumbering]
                            U[i,j]+=localDisp[0]
                            V[i,j]+=localDisp[1]
                            Sxx[i,j]+=localStress[0]
                            Syy[i,j]+=localStress[1]
                            Sxy[i,j]+=localStress[2]
                        elif len(tempCoord)==9:
                            newIsoCoord=invShapeFunction.invQuadQuad(tempCoord,xy)
                            Nmat,Bmat,_=shapeFunction.shapeQuadQuadFE(tempCoord,newIsoCoord[0],newIsoCoord[1])
                            localDisp=rhoValue*Nmat@displacement[fNumbering]
                            localStress=rhoValue*DMatrix@Bmat@displacement[fNumbering]+DMatrix@LCmat@Nmat@displacement[fNumbering]
                            U[i,j]+=localDisp[0]
                            V[i,j]+=localDisp[1]
                            Sxx[i,j]+=localStress[0]
                            Syy[i,j]+=localStress[1]
                            Sxy[i,j]+=localStress[2]
                        else: raise ValueError
        
    elif len(coord)==6:
        newCoord=np.zeros((9,2))
        newCoord[0:3,:]=coord[0:3,:]
        newCoord[6,:]=newCoord[3,:]=coord[2,:]
        newCoord[4,:]=coord[3,:]
        newCoord[5,:]=coord[4,:]
        newCoord[7,:]=newCoord[8,:]=coord[5,:]

        for i in range(nSampling):
            for j in range(nSampling):
                Nmat,_=invShapeFunction.quadQuadShapeFunction([r[i],s[j]])
                xy=(Nmat@newCoord).reshape(-1)
                X[i,j]=xy[0]
                Y[i,j]=xy[1]
                isoCoord=invShapeFunction.invTriQuad(coord,xy)
                triNmat,_=invShapeFunction.triQuadShapeFunction(isoCoord)

                for k in range(len(incidentElements)):
                    if incidentElements[k]:
                        numbering=incidentElements[k].numbering
                        fNumbering=getFullNumber(numbering)
                        tempCoord=AMORE.getCoord(coordinates,numbering)
                        rho=np.zeros(6)
                        rho[0:3]=AMORE.getRho(edge,k)
                        LCmat=stiffnessMatrix.getLC(coord,rho[0:3],isoCoord[0],isoCoord[1])
                        temp=[1,2,0]
                        for l in range(3):
                            rho[l+3]=(rho[l]+rho[temp[l]])/2
                        rhoValue=(triNmat@rho)[0]

                        if len(tempCoord)==4:
                            newIsoCoord=invShapeFunction.invQuad(tempCoord,xy)
                            Nmat,Bmat,_=shapeFunction.shapeQuadFE(tempCoord,newIsoCoord[0],newIsoCoord[1])
                            localDisp=rhoValue*Nmat@displacement[fNumbering]
                            localStress=rhoValue*DMatrix@Bmat@displacement[fNumbering]+DMatrix@LCmat@Nmat@displacement[fNumbering]
                            U[i,j]+=localDisp[0]
                            V[i,j]+=localDisp[1]
                            Sxx[i,j]+=localStress[0]
                            Syy[i,j]+=localStress[1]
                            Sxy[i,j]+=localStress[2]
                        elif len(tempCoord)==9:
                            newIsoCoord=invShapeFunction.invQuadQuad(tempCoord,xy)
                            Nmat,Bmat,_=shapeFunction.shapeQuadQuadFE(tempCoord,newIsoCoord[0],newIsoCoord[1])
                            localDisp=rhoValue*Nmat@displacement[fNumbering]
                            localStress=rhoValue*DMatrix@Bmat@displacement[fNumbering]+DMatrix@LCmat@Nmat@displacement[fNumbering]
                            U[i,j]+=localDisp[0]
                            V[i,j]+=localDisp[1]
                            Sxx[i,j]+=localStress[0]
                            Syy[i,j]+=localStress[1]
                            Sxy[i,j]+=localStress[2]
                        else: raise ValueError

    else: raise ValueError

    return X,Y,U,V,Sxx,Syy,Sxy

def getFullNumber(numbering):
    fNumbering=np.zeros(2*len(numbering),dtype=int)

    for i in range(len(numbering)):
        fNumbering[2*i]=numbering[i]*2
        fNumbering[2*i+1]=numbering[i]*2+1

    return fNumbering

def getGraphDataFE(coord,numbering,displacement,DMatrix,nSampling,solver=3):
    X=np.zeros((nSampling,nSampling))
    Y=np.zeros((nSampling,nSampling))
    U=np.zeros((nSampling,nSampling))
    V=np.zeros((nSampling,nSampling))
    Sxx=np.zeros((nSampling,nSampling))
    Syy=np.zeros((nSampling,nSampling))
    Sxy=np.zeros((nSampling,nSampling))

    r=np.linspace(-1,1,nSampling)
    s=np.linspace(-1,1,nSampling)

    fNumbering=getFullNumber(numbering)

    if len(coord)==3: # Linear triangular element.
        pass
    elif len(coord)==4: # Bilinear quadrilateral element.
        for i in range(nSampling):
            for j in range(nSampling):
                Nmat,_=invShapeFunction.quadShapeFunction([r[i],s[j]])
                xy=Nmat@coord
                Nmat,Bmat,_=shapeFunction.shapeQuadFE(coord,r[i],s[j])
                localDisp=Nmat@displacement[fNumbering]
                if solver==2:
                    localStress=DMatrix@Bmat@displacement[fNumbering]
                elif solver==3: # ICM elements.
                    Bmat,Gmat,_=shapeFunction.shapeICMFE(coord,r[i],s[j])
                    _,Hloc,Eloc=stiffnessMatrix.ICMFE(coord,DMatrix)
                    localStress=DMatrix@(Bmat-Gmat@np.linalg.solve(Hloc,Eloc.transpose()))@displacement[fNumbering]
                else: raise ValueError
                X[i,j]=xy[0,0]
                Y[i,j]=xy[0,1]
                U[i,j]=localDisp[0]
                V[i,j]=localDisp[1]
                Sxx[i,j]=localStress[0]
                Syy[i,j]=localStress[1]
                Sxy[i,j]=localStress[2]

    elif len(coord)==6: # Quadratic triangular element.
        pass
    elif len(coord)==9: # Quadratic quadrilateral element.
        for i in range(nSampling):
            for j in range(nSampling):
                Nmat,_=invShapeFunction.quadQuadShapeFunction([r[i],s[j]])
                xy=Nmat@coord
                Nmat,Bmat,_=shapeFunction.shapeQuadQuadFE(coord,r[i],s[j])
                localDisp=Nmat@displacement[fNumbering]
                localStress=DMatrix@Bmat@displacement[fNumbering]
                X[i,j]=xy[0,0]
                Y[i,j]=xy[0,1]
                U[i,j]=localDisp[0]
                V[i,j]=localDisp[1]
                Sxx[i,j]=localStress[0]
                Syy[i,j]=localStress[1]
                Sxy[i,j]=localStress[2]

    else: raise ValueError

    return X,Y,U,V,Sxx,Syy,Sxy

def getCurvedBoundary(edge,coordinates):
    coord=[]
    startingEdge=edge

    addCoord(coord,edge,coordinates)
    edge=edge.next

    while edge is not startingEdge:
        addCoord(coord,edge,coordinates)
        edge=edge.next    

    return np.array(coord)    

def addCoord(coord,edge,coordinates):
    twinEdge=edge.twin
    numbering=twinEdge.incidentFace.element.numbering

    if len(numbering)>4:
        if len(numbering)==6:
            temp=[1,2,0]
            ind=0
            for i in range(3):
                if AMORE.isEqual(twinEdge.origin,coordinates[numbering[i]]) and \
                    AMORE.isEqual(twinEdge.destination,coordinates[numbering[temp[i]]]):
                    if coordinates[numbering[i+3]] is not None:
                        tempCoord=np.array([coordinates[numbering[temp[i]]],coordinates[numbering[i]],coordinates[numbering[i+3]]])
                        tempCoord=pM.quadInterpolation(tempCoord,5)
                        for j in range(1,len(tempCoord)):
                            coord.append([tempCoord[j,0],tempCoord[j,1]])
                    else: coord.append([edge.destination.x,edge.destination.y])
                    ind=1
                    break

            if ind==0: coord.append([edge.destination.x,edge.destination.y])
        elif len(numbering)==9:
            temp=[1,2,3,0]
            ind=0
            for i in range(4):
                if AMORE.isEqual(twinEdge.origin,coordinates[numbering[i]]) and \
                    AMORE.isEqual(twinEdge.destination,coordinates[numbering[temp[i]]]):
                    if coordinates[numbering[i+4]] is not None:
                        tempCoord=np.array([coordinates[numbering[temp[i]]],coordinates[numbering[i]],coordinates[numbering[i+4]]])
                        tempCoord=pM.quadInterpolation(tempCoord,5)
                        for j in range(1,len(tempCoord)):
                            coord.append([tempCoord[j,0],tempCoord[j,1]])
                    else: coord.append([edge.destination.x,edge.destination.y])
                    ind=1
                    break

            if ind==0: coord.append([edge.destination.x,edge.destination.y])
        else: raise ValueError

    else: coord.append([edge.destination.x,edge.destination.y]) 

def getNumbering(incidentElements):
    for i in incidentElements:
        if i: return i.numbering
    
    return None

if __name__=='__main__':
    pass
