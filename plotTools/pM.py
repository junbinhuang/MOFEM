import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np
import math, copy
from matplotlib import rc
import sys, os
sys.path.append(os.path.dirname(sys.path[0]))
from linearSolvers import AMORE
from elementLibrary import shapeFunction
from plotTools import plotResults

######
# Function List:
# plotMesh
######

def plotMesh(coordinates,meshes,fmt,lineWidth,markerSize):
    """Plot a two-dimensional mesh."""
    # fmt: The format string
    # lineWidth
    # markerSize

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    if len(meshes)==1:
        elements=meshes[0]
        for i in elements:
            coord=AMORE.getCoord(coordinates,i)
            plotElement(coord,fmt,lineWidth,markerSize)
    else:
        colorList=['b','r','g','c','m']
        assert len(colorList)>=len(meshes) # I don't see any need for using more than 5 meshes.

        count=0
        for i in range(len(colorList)):
            if colorList[i]==fmt[0]: 
                count=i
                break
        if count!=0: colorList[0],colorList[count]=colorList[count],colorList[0]

        fmtList=[None]*len(meshes)
        fmtList[0]=fmt
        for i in range(1,len(fmtList)):
            fmtList[i]=colorList[i]+fmt[1:len(fmt)]

        for k in range(len(meshes)):
            newMarkerSize=markerSize*(0.7)**(k>0)
            elements=meshes[k]

            for i in elements:
                coord=AMORE.getCoord(coordinates,i)
                plotElement(coord,fmtList[k],lineWidth,newMarkerSize)

    plt.axis('off')
    plt.axis('equal')
    plt.title('Input Mesh')
    plt.savefig('Input Mesh.pdf',bbox_inches='tight')
    plt.show()

def findCoord(number,coordinates,meshes):
    x=y=None

    for mesh in meshes:
        for i in mesh:
            for j in range(len(i)):
                if i[j]==number: 
                    if len(i)==6:
                        array=[1,2,0]
                        x=(coordinates[i[j-3]][0]+coordinates[i[array[j-3]]][0])/2
                        y=(coordinates[i[j-3]][1]+coordinates[i[array[j-3]]][1])/2
                    elif len(i)==9:
                        assert(j<8),"It is a central node!"
                        array=[1,2,3,0]
                        x=(coordinates[i[j-4]][0]+coordinates[i[array[j-4]]][0])/2
                        y=(coordinates[i[j-4]][1]+coordinates[i[array[j-4]]][1])/2
                    else: raise ValueError

    assert(x is not None and y is not None),"Cannot find this node!"

    return x,y

def plotElement(coord,fmt,lineWidth,markerSize):
    if len(coord)==3 or len(coord)==4: # 1st-order element.
        plt.plot(coord[:,0],coord[:,1],fmt,linewidth=lineWidth,markersize=markerSize,fillstyle='full')
        plt.plot(coord[[0,-1],0],coord[[0,-1],1],fmt,linewidth=lineWidth,markersize=markerSize,fillstyle='full')
    elif len(coord)==6: # 2nd-order triangle.
        fmt1=fmt[0]+fmt[2]
        fmt2=fmt[0]+fmt[1]
        temp=[1,2,0]

        for i in range(3):
            newCoord=quadInterpolation(coord[[i,temp[i],i+3],:],7)
            plt.plot(newCoord[:,0],newCoord[:,1],fmt1,linewidth=lineWidth,markersize=markerSize,fillstyle='full')
        plt.plot(coord[:,0],coord[:,1],fmt2,linewidth=lineWidth,markersize=markerSize,fillstyle='full')
    elif len(coord)==9: # 2nd-order quad.
        fmt1=fmt[0]+fmt[2]
        fmt2=fmt[0]+fmt[1]
        temp=[1,2,3,0]

        for i in range(4):
            newCoord=quadInterpolation(coord[[i,temp[i],i+4],:],7)
            plt.plot(newCoord[:,0],newCoord[:,1],fmt1,linewidth=lineWidth,markersize=markerSize,fillstyle='full')
        plt.plot(coord[:,0],coord[:,1],fmt2,linewidth=lineWidth,markersize=markerSize,fillstyle='full')

def quadInterpolation(coord,n):
    """Input:
    
    coord=array([[x1,y1],[x2,y2],[x3,y3]]), where [x3,y3] is the mid-point."""
    newCoord=np.zeros((n,2))
    for i in range(n):
        r=-1.0+2.0*i/(n-1.0)
        Nmat=shapeFunction.oneDQuadratic2(r)
        newCoord[i]=(Nmat@coord).reshape(-1)
    
    return newCoord