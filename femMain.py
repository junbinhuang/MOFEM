######
# FEM in Python
######

######
# Import modules
import numpy as np
from meshTools import toDC
from inputFunctions import inputF
from outputFunctions import writeSolution
from plotTools import pM, plotResults
from computGeometry import planeSweepMeshes, triangulation
from linearSolvers import traditionalElement, AMORE
import time, copy
######

######
if __name__=="__main__":
    # Read data from a typical FE input file
    inputFileName=input("Please input the file name: ")
    # Warning: For large problems, give coodinates small perturbations before running the plane sweep algorithm.

    startTime=time.time()
    inputData=inputF.txtInput(inputFileName)
    meshes=inputData[3]
    coordinates=inputData[5]
    constraints=inputData[7]
    print("Importing data costs %s seconds."%(time.time()-startTime))

    # pM.plotMesh(coordinates,meshes,"bo-",1,2) # Plot meshes

    startTime=time.time()
    parameters=inputData[0]
    numberOfMesh=len(meshes)
    InputMeshes=copy.deepcopy(meshes)
    _,tempMesh=planeSweepMeshes.elementFilter(coordinates,InputMeshes,transform=True)
    # Warning: For large problems, give coodinates small perturbations before running the plane sweep algorithm.
    tempVertexList,tempFaceList,tempHalfEdgeList,tempLineList,tempNodeElementsPartial=toDC.toDoublyConnected(coordinates,tempMesh)
    # This is used for visualization only, see the Matlab code.
    tempPolygon,boundary=planeSweepMeshes.meshOverlay(tempLineList,tempVertexList,tempHalfEdgeList,tempFaceList,numberOfMesh)

    for i in tempPolygon:
        triangulation.polygonTriangulation(i)

    print("Getting boundary costs %s seconds."%(time.time()-startTime))

    startTime=time.time()

    if len(meshes)==1:
        # Only one mesh
        polygons=[]
        if parameters[0]==1:
            # 3-node triangular FE or 4-node quadrilateral FE
            print("Using traditional lower-order finite elements.")
            displacement,energy=traditionalElement.lowerOrderFE(inputData)
            print("FE solver costs %s seconds in total."%(time.time()-startTime))
            writeSolution.writeDisplacement(displacement,energy,inputFileName)

        elif parameters[0]==2:
            # 6-node triangular FE or 9-node quadrilateral FE
            print("Using traditional quadratic finite elements.")
            displacement,energy=traditionalElement.quadraticFE(inputData)
            print("FE solver costs %s seconds in total."%(time.time()-startTime))
            writeSolution.writeDisplacement(displacement,energy,inputFileName)

        elif parameters[0]==3:
            # ICM FE (for squares)
            print("Using incompatible finite elements.")
            displacement,energy=traditionalElement.ICMFE(inputData)
            print("FE solver costs %s seconds in total."%(time.time()-startTime))
            writeSolution.writeDisplacement(displacement,energy,inputFileName)

        else:
            raise ValueError("Wrong element type!")
    else:
        # Handle the overlapping meshes.
        (meshes,overlappingMesh)=planeSweepMeshes.elementFilter(coordinates,meshes)
        vertexList,faceList,halfEdgeList,lineList,nodeElementsPartial=toDC.toDoublyConnected(coordinates,overlappingMesh)
        # Warning: For large problems, give coodinates small perturbations before running the plane sweep algorithm.
        polygons,_=planeSweepMeshes.meshOverlay(lineList,vertexList,halfEdgeList,faceList,numberOfMesh)
        planeSweepMeshes.rhoCalculation(polygons)

        for i in polygons:
            triangulation.polygonTriangulation(i)
                
        print("Overlapping calculation costs %s seconds."%(time.time()-startTime))
        startTime=time.time()

        if parameters[0]==1:
            # We use 3-node FE or 4-node FE for the boundary meshes.
            print("Using 4-node finite elements for the interior and lower-order finite elements for the boundary meshes.")
            displacement,energy=AMORE.lowerOrderAMORE(inputData,overlappingMesh,polygons,nodeElementsPartial)
            print("FE solver costs %s seconds in total."%(time.time()-startTime))
            writeSolution.writeDisplacement(displacement,energy,inputFileName)

        elif parameters[0]==2:
            # We use 6-node FE or 9-node FE for the boundary meshes.
            print("Using 4-node finite elements for the interior and quadratic finite elements for the boundary meshes.")
            displacement,energy,materialMatrix=AMORE.quadraticAMORE(inputData,overlappingMesh,polygons,nodeElementsPartial)
            print("FE solver costs %s seconds in total."%(time.time()-startTime))
            writeSolution.writeDisplacement(displacement,energy,inputFileName)
            writeSolution.writeBoundary(boundary,coordinates)
            writeSolution.writeSampleData2D(coordinates,meshes,polygons,displacement,materialMatrix,inputData[4],component='All',nSampling=3,solver=2)

        elif parameters[0]==3:
            # We use 6-node FE or 9-node FE for the boundary meshes. And ICM element is used to replace all regular, non-overlapping 4-node elements.
            print("Using ICM elements for the interior and quadratic finite elements for the boundary meshes.")
            displacement,energy,materialMatrix=AMORE.quadraticICMAMORE(inputData,overlappingMesh,polygons,nodeElementsPartial)
            print("FE solver costs %s seconds in total."%(time.time()-startTime))
            writeSolution.writeDisplacement(displacement,energy,inputFileName)
            writeSolution.writeBoundary(boundary,coordinates)
            writeSolution.writeSampleData2D(coordinates,meshes,polygons,displacement,materialMatrix,inputData[4],component='All',nSampling=3,solver=3)
        else:
            raise ValueError("Wrong element type!")