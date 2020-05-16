import sys, os
import numpy as np
sys.path.append(os.path.dirname(sys.path[0]))
from plotTools import plotResults
from linearSolvers import AMORE
import math

def writeDisplacement(displacement,energy=None,inputFileName=None):
    file=open("solution.txt","w")
    if inputFileName is not None:
        file.write("The input file name is %s.txt.\n\n"%inputFileName)

    file.write("                      DISPLACEMENT\n\n")
    file.write("  Node No.           x-displacement                   y-displacement\n")
    for i in range(displacement.size//2):
        file.write("%6d           %22.15e           %22.15e\n"%(i,displacement[2*i][0],displacement[2*i+1][0]))
    
    file.write("\n")
    if energy is not None: 
        file.write("The strain ENERGY is: %.15e."%energy[0][0])

def writeRow(list,file):
    for i in list: file.write("%s "%i)
    file.write("\n")

def writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file):
    for k1 in range(nSampling):
        writeRow(X[k1],file)
        writeRow(Y[k1],file)
        writeRow(Z[k1],file)
        for k2 in range(nSampling):
            if Z[k1,k2]>limit[1]: 
                limit[1]=Z[k1,k2]
                maxCoord[0]=X[k1,k2]
                maxCoord[1]=Y[k1,k2]
            if Z[k1,k2]<limit[0]: 
                limit[0]=Z[k1,k2]
                minCoord[0]=X[k1,k2]
                minCoord[1]=Y[k1,k2]

def writeSampleData2D(coordinates,meshes,polygons,displacement,materialMatrix,materialMeshList,component='u',nSampling=None,solver=3):
    """If component='u', write down the u-displacements;
       If component='v', write down the v-displacements;
       If component='Sxx', write down stress_xx;
       ...          'Sxy', ...        stress_xy;
       ...          'Syy', ...        stress_yy;
       ...          'All', ...        all displacement and stress components.
       """
    if nSampling is None: nSampling=5
    file=open("SampleData.txt","w")

    print("%s sampling points are used in each direction."%nSampling)
    limit=[math.inf,-math.inf]
    minCoord=[0.0,0.0]
    maxCoord=[0.0,0.0]

    if component=='u':
        for i in range(len(meshes)):
            for j in range(len(meshes[i])):
                if meshes[i][j]:
                    DMatrix=materialMatrix[materialMeshList[i][j]]
                    coord=AMORE.getCoord(coordinates,meshes[i][j])
                    X,Y,Z,_,_,_,_=plotResults.getGraphDataFE(coord,meshes[i][j],displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

        if polygons:
            for i in polygons:
                if i.triangles:
                    for j in i.triangles:
                        for k in i.incidentElements:
                            if k:
                                position=k.position
                                break
                        
                        DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                        coord=AMORE.coordCurvedTriangle(j,i.incidentElements,coordinates)
                        X,Y,Z,_,_,_,_=plotResults.getGraphDataOFE(coord,coordinates,i.incidentElements,j,displacement,DMatrix,nSampling)
                        writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)
                else:
                    for k in i.incidentElements:
                        if k:
                            position=k.position
                            break

                    DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                    numbering=plotResults.getNumbering(i.incidentElements)
                    coord=AMORE.getCoord(coordinates,numbering)
                    X,Y,Z,_,_,_,_=plotResults.getGraphDataFE(coord,numbering,displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

    elif component=='v':
        for i in range(len(meshes)):
            for j in range(len(meshes[i])):
                if meshes[i][j]:
                    DMatrix=materialMatrix[materialMeshList[i][j]]
                    coord=AMORE.getCoord(coordinates,meshes[i][j])
                    X,Y,_,Z,_,_,_=plotResults.getGraphDataFE(coord,meshes[i][j],displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

        if polygons:
            for i in polygons:
                if i.triangles:
                    for j in i.triangles:
                        for k in i.incidentElements:
                            if k:
                                position=k.position
                                break
                        
                        DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                        coord=AMORE.coordCurvedTriangle(j,i.incidentElements,coordinates)
                        X,Y,_,Z,_,_,_=plotResults.getGraphDataOFE(coord,coordinates,i.incidentElements,j,displacement,DMatrix,nSampling)
                        writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)
                else:
                    for k in i.incidentElements:
                        if k:
                            position=k.position
                            break

                    DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]
                    
                    numbering=plotResults.getNumbering(i.incidentElements)
                    coord=AMORE.getCoord(coordinates,numbering)
                    X,Y,_,Z,_,_,_=plotResults.getGraphDataFE(coord,numbering,displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)
    
    elif component=='Sxx':
        for i in range(len(meshes)):
            for j in range(len(meshes[i])):
                if meshes[i][j]:
                    DMatrix=materialMatrix[materialMeshList[i][j]]
                    coord=AMORE.getCoord(coordinates,meshes[i][j])
                    X,Y,_,_,Z,_,_=plotResults.getGraphDataFE(coord,meshes[i][j],displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

        if polygons:
            for i in polygons:
                if i.triangles:
                    for j in i.triangles:
                        for k in i.incidentElements:
                            if k:
                                position=k.position
                                break
                        
                        DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                        coord=AMORE.coordCurvedTriangle(j,i.incidentElements,coordinates)
                        X,Y,_,_,Z,_,_=plotResults.getGraphDataOFE(coord,coordinates,i.incidentElements,j,displacement,DMatrix,nSampling)
                        writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)
                else:
                    for k in i.incidentElements:
                        if k:
                            position=k.position
                            break

                    DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]
                    
                    numbering=plotResults.getNumbering(i.incidentElements)
                    coord=AMORE.getCoord(coordinates,numbering)
                    X,Y,_,_,Z,_,_=plotResults.getGraphDataFE(coord,numbering,displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

    elif component=='Sxy':
        for i in range(len(meshes)):
            for j in range(len(meshes[i])):
                if meshes[i][j]:
                    DMatrix=materialMatrix[materialMeshList[i][j]]
                    coord=AMORE.getCoord(coordinates,meshes[i][j])
                    X,Y,_,_,_,_,Z=plotResults.getGraphDataFE(coord,meshes[i][j],displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

        if polygons:
            for i in polygons:
                if i.triangles:
                    for j in i.triangles:
                        for k in i.incidentElements:
                            if k:
                                position=k.position
                                break
                        
                        DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                        coord=AMORE.coordCurvedTriangle(j,i.incidentElements,coordinates)
                        X,Y,_,_,_,_,Z=plotResults.getGraphDataOFE(coord,coordinates,i.incidentElements,j,displacement,DMatrix,nSampling)
                        writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)
                else:
                    for k in i.incidentElements:
                        if k:
                            position=k.position
                            break

                    DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]
                    
                    numbering=plotResults.getNumbering(i.incidentElements)
                    coord=AMORE.getCoord(coordinates,numbering)
                    X,Y,_,_,_,_,Z=plotResults.getGraphDataFE(coord,numbering,displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

    elif component=='Syy':
        for i in range(len(meshes)):
            for j in range(len(meshes[i])):
                if meshes[i][j]:
                    DMatrix=materialMatrix[materialMeshList[i][j]]
                    coord=AMORE.getCoord(coordinates,meshes[i][j])
                    X,Y,_,_,_,Z,_=plotResults.getGraphDataFE(coord,meshes[i][j],displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

        if polygons:
            for i in polygons:
                if i.triangles:
                    for j in i.triangles:
                        for k in i.incidentElements:
                            if k:
                                position=k.position
                                break
                        
                        DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                        coord=AMORE.coordCurvedTriangle(j,i.incidentElements,coordinates)
                        X,Y,_,_,_,Z,_=plotResults.getGraphDataOFE(coord,coordinates,i.incidentElements,j,displacement,DMatrix,nSampling)
                        writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)
                else:
                    for k in i.incidentElements:
                        if k:
                            position=k.position
                            break

                    DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]
                    
                    numbering=plotResults.getNumbering(i.incidentElements)
                    coord=AMORE.getCoord(coordinates,numbering)
                    X,Y,_,_,_,Z,_=plotResults.getGraphDataFE(coord,numbering,displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limit,minCoord,maxCoord,X,Y,Z,nSampling,file)

    elif component=='All':
        limitU=[math.inf,-math.inf]
        minCoordU=[0.0,0.0]
        maxCoordU=[0.0,0.0]
        limitV=[math.inf,-math.inf]
        minCoordV=[0.0,0.0]
        maxCoordV=[0.0,0.0]
        limitSxx=[math.inf,-math.inf]
        minCoordSxx=[0.0,0.0]
        maxCoordSxx=[0.0,0.0]
        limitSyy=[math.inf,-math.inf]
        minCoordSyy=[0.0,0.0]
        maxCoordSyy=[0.0,0.0]
        limitSxy=[math.inf,-math.inf]
        minCoordSxy=[0.0,0.0]
        maxCoordSxy=[0.0,0.0]
        limitMises=[math.inf,-math.inf]
        minCoordMises=[0.0,0.0]
        maxCoordMises=[0.0,0.0]

        for i in range(len(meshes)):
            for j in range(len(meshes[i])):
                if meshes[i][j]:
                    DMatrix=materialMatrix[materialMeshList[i][j]]
                    coord=AMORE.getCoord(coordinates,meshes[i][j])
                    X,Y,U,V,Sxx,Syy,Sxy=plotResults.getGraphDataFE(coord,meshes[i][j],displacement,DMatrix,nSampling,solver)

                    writeAndUpdate(limitU,minCoordU,maxCoordU,X,Y,U,nSampling,file)
                    writeAndUpdate(limitV,minCoordV,maxCoordV,X,Y,V,nSampling,file)
                    writeAndUpdate(limitSxx,minCoordSxx,maxCoordSxx,X,Y,Sxx,nSampling,file)
                    writeAndUpdate(limitSyy,minCoordSyy,maxCoordSyy,X,Y,Syy,nSampling,file)
                    writeAndUpdate(limitSxy,minCoordSxy,maxCoordSxy,X,Y,Sxy,nSampling,file)

                    Mises=np.zeros((nSampling,nSampling))
                    for ii in range(nSampling):
                        for jj in range(nSampling):
                            Mises[ii,jj]=math.sqrt(Sxx[ii,jj]**2+Syy[ii,jj]**2+3*Sxy[ii,jj]**2-Sxx[ii,jj]*Syy[ii,jj])
                    writeAndUpdate(limitMises,minCoordMises,maxCoordMises,X,Y,Mises,nSampling,file)

        if polygons:
            for i in polygons:
                if i.triangles:
                    for j in i.triangles:
                        for k in i.incidentElements:
                            if k:
                                position=k.position
                                break
                        
                        DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                        coord=AMORE.coordCurvedTriangle(j,i.incidentElements,coordinates)
                        if abs(triArea(coord))<1e-6: continue # These small triangles may bring much error to the stress results.

                        X,Y,U,V,Sxx,Syy,Sxy=plotResults.getGraphDataOFE(coord,coordinates,i.incidentElements,j,displacement,DMatrix,nSampling)
                        writeAndUpdate(limitU,minCoordU,maxCoordU,X,Y,U,nSampling,file)
                        writeAndUpdate(limitV,minCoordV,maxCoordV,X,Y,V,nSampling,file)
                        writeAndUpdate(limitSxx,minCoordSxx,maxCoordSxx,X,Y,Sxx,nSampling,file)
                        writeAndUpdate(limitSyy,minCoordSyy,maxCoordSyy,X,Y,Syy,nSampling,file)
                        writeAndUpdate(limitSxy,minCoordSxy,maxCoordSxy,X,Y,Sxy,nSampling,file)

                        Mises=np.zeros((nSampling,nSampling))
                        for ii in range(nSampling):
                            for jj in range(nSampling):
                                Mises[ii,jj]=math.sqrt(Sxx[ii,jj]**2+Syy[ii,jj]**2+3*Sxy[ii,jj]**2-Sxx[ii,jj]*Syy[ii,jj])
                        writeAndUpdate(limitMises,minCoordMises,maxCoordMises,X,Y,Mises,nSampling,file)

                else:
                    for k in i.incidentElements:
                        if k:
                            position=k.position
                            break
                    
                    DMatrix=materialMatrix[materialMeshList[position[0]][position[1]]]

                    numbering=plotResults.getNumbering(i.incidentElements)
                    coord=AMORE.getCoord(coordinates,numbering)
                    X,Y,U,V,Sxx,Syy,Sxy=plotResults.getGraphDataFE(coord,numbering,displacement,DMatrix,nSampling,solver)
                    writeAndUpdate(limitU,minCoordU,maxCoordU,X,Y,U,nSampling,file)
                    writeAndUpdate(limitV,minCoordV,maxCoordV,X,Y,V,nSampling,file)
                    writeAndUpdate(limitSxx,minCoordSxx,maxCoordSxx,X,Y,Sxx,nSampling,file)
                    writeAndUpdate(limitSyy,minCoordSyy,maxCoordSyy,X,Y,Syy,nSampling,file)
                    writeAndUpdate(limitSxy,minCoordSxy,maxCoordSxy,X,Y,Sxy,nSampling,file)

                    Mises=np.zeros((nSampling,nSampling))
                    for ii in range(nSampling):
                        for jj in range(nSampling):
                            Mises[ii,jj]=math.sqrt(Sxx[ii,jj]**2+Syy[ii,jj]**2+3*Sxy[ii,jj]**2-Sxx[ii,jj]*Syy[ii,jj])
                    writeAndUpdate(limitMises,minCoordMises,maxCoordMises,X,Y,Mises,nSampling,file)

    else: raise ValueError
    if component!='All':
        print("The range of "+component+" is [%s,%s]."%(limit[0],limit[1]))
        print("Location of the minimum is [%s,%s]."%(minCoord[0],minCoord[1]))
        print("Location of the maximum is [%s,%s]."%(maxCoord[0],maxCoord[1]))
    else:
        print("The range of u is [%s,%s]."%(limitU[0],limitU[1]))
        print("Location of the minimum is [%s,%s]."%(minCoordU[0],minCoordU[1]))
        print("Location of the maximum is [%s,%s]."%(maxCoordU[0],maxCoordU[1]))
        print("The range of v is [%s,%s]."%(limitV[0],limitV[1]))
        print("Location of the minimum is [%s,%s]."%(minCoordV[0],minCoordV[1]))
        print("Location of the maximum is [%s,%s]."%(maxCoordV[0],maxCoordV[1]))
        print("The range of Sxx is [%s,%s]."%(limitSxx[0],limitSxx[1]))
        print("Location of the minimum is [%s,%s]."%(minCoordSxx[0],minCoordSxx[1]))
        print("Location of the maximum is [%s,%s]."%(maxCoordSxx[0],maxCoordSxx[1]))
        print("The range of Syy is [%s,%s]."%(limitSyy[0],limitSyy[1]))
        print("Location of the minimum is [%s,%s]."%(minCoordSyy[0],minCoordSyy[1]))
        print("Location of the maximum is [%s,%s]."%(maxCoordSyy[0],maxCoordSyy[1]))
        print("The range of Sxy is [%s,%s]."%(limitSxy[0],limitSxy[1]))
        print("Location of the minimum is [%s,%s]."%(minCoordSxy[0],minCoordSxy[1]))
        print("Location of the maximum is [%s,%s]."%(maxCoordSxy[0],maxCoordSxy[1]))
        print("The range of Mises stress is [%s,%s]."%(limitMises[0],limitMises[1]))
        print("Location of the minimum is [%s,%s]."%(minCoordMises[0],minCoordMises[1]))
        print("Location of the maximum is [%s,%s]."%(maxCoordMises[0],maxCoordMises[1]))        

    file=open("nSampling.txt","w")
    file.write(str(nSampling))
    file.write("\n")
    if component!='All': file.write(str(1))
    else: file.write(str(6)) # On Dec 11th, changed 5 to 6 since we want to output Mises stress as well.

def writeBoundary(edgeList,coordinates):
    length=[]
    file=open("boundaryCoord.txt","w")

    for edge in edgeList:
        coord=plotResults.getCurvedBoundary(edge,coordinates)
        for i in coord:
            writeRow(i,file)
        length.append(len(coord))

    file=open("boundaryNumber.txt","w")
    for i in length:
        file.write("%s\n"%i)

def triArea(coord):
    """Input:
    
    coord=array([[x1,y1],[x2,y2],[x3,y3]])."""
    area=(coord[1,0]*coord[2,1]-coord[2,0]*coord[1,1]-coord[0,0]*\
        (coord[2,1]-coord[1,1])+coord[0,1]*(coord[2,0]-coord[1,0]))/2
    
    return area

if __name__=="__main__":
    pass