# The AMORE scheme for linear problems.
import numpy as np 
import math
import scipy
# from sksparse.cholmod import cholesky # It works.
from scipy.sparse.linalg import spsolve
import time
import sys, os
sys.path.append(os.path.dirname(sys.path[0]))
from elementLibrary import stiffnessMatrix, shapeFunction, invShapeFunction
from otherFunctions import numericalIntegration
from computGeometry import planeSweepMeshes, triangulation
from meshTools import toDC
from linearSolvers import traditionalElement
import cForce

def lowerOrderAMORE(inputData,overlappingMesh,polygons,nodeElementsPartial):
    """AMORE scheme with lower-order traditional FE on the boundary meshes."""
    startTime=time.time()

    parameters=inputData[0]

    # Material matrices.
    materialList=inputData[6]

    materialMatrix=[None]*len(materialList)

    for i in range(len(materialList)):
        materialMatrix[i]=twoDMaterialMatrix(materialList[i],parameters[1])
    
    # Assemble stiffness matrix.
    coordinates=inputData[5]
    meshes=inputData[3]
    numberOfMesh=len(meshes)
    materialMeshList=inputData[4]

    Iglo=[]
    Jglo=[]
    Vglo=[]

    # Integrate over non-overlapping elements.
    for i in range(len(meshes)):
        for j in range(len(meshes[i])):
            if meshes[i][j] is not None:
                coord=getCoord(coordinates,meshes[i][j])    

                if len(meshes[i][j])==3:
                    Kloc=stiffnessMatrix.triFE(coord,materialMatrix[materialMeshList[i][j]])
                elif len(meshes[i][j])==4:
                    Kloc=stiffnessMatrix.quadFE(coord,materialMatrix[materialMeshList[i][j]])
                else: raise ValueError("Wrong element numbering!")
            
                Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,meshes[i][j])
                Iglo.extend(Iloc)
                Jglo.extend(Jloc)
                Vglo.extend(Vloc)
    
    # Integrate for the overlapping meshes.
    for i in range(len(polygons)):
        if polygons[i].triangles is None:
            # It is a traditional element.
            for j in polygons[i].incidentElements:
                if j:
                    numbering=j.numbering
                    position=j.position
                    break
            
            coord=getCoord(coordinates,numbering)
            
            if len(numbering)==3:
                Kloc=stiffnessMatrix.triFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            elif len(numbering)==4:
                Kloc=stiffnessMatrix.quadFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            else: raise ValueError("Wrong element numbering!")

            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering)
            Iglo.extend(Iloc)
            Jglo.extend(Jloc)
            Vglo.extend(Vloc)
        else:
            for j in range(numberOfMesh):
                if polygons[i].incidentElements[j]:
                    numbering=polygons[i].incidentElements[j].numbering
                    position=polygons[i].incidentElements[j].position
                    coord=getCoord(coordinates,numbering)
                    Kloc=np.zeros((2*len(numbering),2*len(numbering)))

                    # Integrate the diagonal terms.
                    for k in polygons[i].triangles:
                        rho=getRho(k,position[0])
                        coordTri=stiffnessMatrix.getCoordFromHalfEdge(k)
                        Kloc+=stiffnessMatrix.lowerOrderAMORE(coordTri,\
                            materialMatrix[materialMeshList[position[0]][position[1]]],coord,rho)

                    Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering)
                    Iglo.extend(Iloc)
                    Jglo.extend(Jloc)
                    Vglo.extend(Vloc)
                        
                    # Integrate the off-diagonal terms.
                    for l in range(j+1,numberOfMesh):
                        if polygons[i].incidentElements[l]:
                            numbering2=polygons[i].incidentElements[l].numbering
                            position2=polygons[i].incidentElements[l].position
                            coord2=getCoord(coordinates,numbering2)
                            Kloc=np.zeros((2*len(numbering),2*len(numbering2)))

                            for k in polygons[i].triangles:
                                rho=getRho(k,position[0])
                                rho2=getRho(k,position2[0])
                                coordTri=stiffnessMatrix.getCoordFromHalfEdge(k)
                                Kloc+=stiffnessMatrix.lowerOrderAMORE(coordTri,\
                                    materialMatrix[materialMeshList[position[0]][position[1]]],coord,rho,coord2,rho2)

                            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering,numbering2)
                            Iglo.extend(Iloc)
                            Jglo.extend(Jloc)
                            Vglo.extend(Vloc)

                            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc.transpose(),numbering2,numbering)
                            Iglo.extend(Iloc)
                            Jglo.extend(Jloc)
                            Vglo.extend(Vloc)                                

    Iglo=np.array(Iglo,dtype=int)
    Jglo=np.array(Jglo,dtype=int)
    Vglo=np.array(Vglo,dtype='d')

    Kglo=scipy.sparse.coo_matrix((Vglo,(Iglo,Jglo)),shape=(2*len(coordinates),2*len(coordinates))).tocsr()
    
    print("Assembling stiffness matrix costs %s seconds."%(time.time()-startTime))
    startTime=time.time()

    # Force term.
    indForce=inputData[-2]
    if indForce[0]: # Body force is imposed.
        pass

    assert (indForce[1]==0),"Customized force is not allowed in this solver!"

    fglo=np.zeros((2*len(coordinates),1))

    if len(indForce)==2 or (len(indForce)==3 and indForce[2]==0):
        forceList=inputData[-1]
    elif len(indForce)==3 and indForce[2]==1:
        forceList=inputData[-1][0]
        cForceList=inputData[-1][1]

        for i in cForceList: # We assume here that no concentrated force is acting on the overlap.
            # We can remove this assumption but the code would be complicated as we need the rho functions.
            fglo[2*i[0],0]+=i[1]
            fglo[2*i[0]+1,0]+=i[2]

    nInt=2
    if indForce[1]: nInt+=3 # Customized boundary force.

    pos,wei=numericalIntegration.gaussQuad(nInt)

    for i in range(len(forceList)//2):
        node1=forceList[2*i][0]
        node2=forceList[2*i+1][0]

        force1=np.array([forceList[2*i][1:3]]).transpose()
        force2=np.array([forceList[2*i+1][1:3]]).transpose()

        floc=np.zeros((4,1))

        coord=np.array([coordinates[node1],coordinates[node2]])
        length=lenEdge(coord[0],coord[1])
        # Find the element.
        elementPosition=None
        element=(set(nodeElementsPartial[node1]) & set(nodeElementsPartial[node2]))
        if element: elementPosition=element.pop()

        if elementPosition and (overlappingMesh[elementPosition[0]]\
            [elementPosition[1]].polygons[0].triangles is not None):
            # It is in some overlapping element and this element is triangulated.
            for j in overlappingMesh[elementPosition[0]][elementPosition[1]].polygons:
                for k in j.triangles:
                    # Test if this triangle is on the boundary.
                    # If yes, integrate.
                    boundaryEdge=isTriangleOnBoundary(k,coord[0],coord[1])
                    if boundaryEdge:
                        rho=np.array([boundaryEdge.origin.rho[elementPosition[0]],\
                            boundaryEdge.destination.rho[elementPosition[0]]])
                        coordEdge=np.array([[boundaryEdge.origin.x,boundaryEdge.origin.y],\
                            [boundaryEdge.destination.x,boundaryEdge.destination.y]])
                        lengthEdge=lenEdge(coordEdge[0],coordEdge[1])

                        for l in range(nInt):
                            Nmat=shapeFunction.oneDLinear2(pos[l])
                            rhoValue=(Nmat@rho)[0]
                            xy=(Nmat@coordEdge).reshape(-1)
                            isoCoord=invShapeFunction.invLine(coord,xy)
                            Nmat=shapeFunction.oneDLinear(isoCoord)
                            force=Nmat[0,0]*force1+Nmat[0,2]*force2
                            floc+=0.5*wei[l]*lengthEdge*rhoValue*np.matmul(Nmat.transpose(),force)
                        break

        else: # This is a regular element.
            for j in range(nInt):
                Nmat=shapeFunction.oneDLinear(pos[j])
                force=Nmat[0,0]*force1+Nmat[0,2]*force2

                floc+=0.5*wei[j]*length*np.matmul(Nmat.transpose(),force)

        fglo[2*node1:2*node1+2,0]+=floc[0:2,0]
        fglo[2*node2:2*node2+2,0]+=floc[2:4,0]
    
    print("Calculating force term costs %s seconds."%(time.time()-startTime))
    startTime=time.time()
        
    # Impose constraints.
    fixList=np.zeros((2*len(coordinates),1))
    fixIndexList=np.zeros((2*len(coordinates),1),dtype=int)

    constraintList=inputData[-3]
    # Very important!!! Sort the constraints!!!
    constraintList.sort(key=lambda item:item[0])

    for i in constraintList:
        if i[1]: 
            fixList[2*i[0]]=i[3]
            fixIndexList[2*i[0]]=1
        if i[2]: 
            fixList[2*i[0]+1]=i[4]
            fixIndexList[2*i[0]+1]=1

    # Solve.
    fglo-=(Kglo.dot(fixList))

    Kglo_complete=Kglo.copy()
    Kglo=Kglo.tolil()

    count=0
    for i in constraintList:
        if i[1]:
            delete_row_lil(Kglo,2*i[0]-count)
            fglo=np.delete(fglo,2*i[0]-count)
            count+=1
        if i[2]:
            delete_row_lil(Kglo,2*i[0]+1-count)
            fglo=np.delete(fglo,2*i[0]+1-count)
            count+=1

    Kglo=Kglo.transpose()

    count=0
    for i in constraintList:
        if i[1]:
            delete_row_lil(Kglo,2*i[0]-count)
            count+=1
        if i[2]:
            delete_row_lil(Kglo,2*i[0]+1-count)
            count+=1

    print("Imposing constraints costs %s seconds."%(time.time()-startTime))
    startTime=time.time()
    
    Kglo=Kglo.tocsc()

    # factor=cholesky(Kglo)
    # disp=factor(fglo)
    disp=spsolve(Kglo,fglo)

    print("Solving the linear system costs %s seconds."%(time.time()-startTime))

    # The complete displacement solution:
    displacement=np.zeros((2*len(coordinates),1))
    count=0
    for i in range(2*len(coordinates)):
        if fixIndexList[i]: 
            displacement[i]=fixList[i]
            count+=1
        else: 
            displacement[i]=disp[i-count]

    energy=0.5*displacement.transpose()@Kglo_complete@displacement

    return displacement,energy

def quadraticAMORE(inputData,overlappingMesh,polygons,nodeElementsPartial):
    """AMORE scheme with lower-order traditional FE on the boundary meshes."""
    startTime=time.time()

    parameters=inputData[0]

    # Material matrices.
    materialList=inputData[6]

    materialMatrix=[None]*len(materialList)

    for i in range(len(materialList)):
        materialMatrix[i]=twoDMaterialMatrix(materialList[i],parameters[1])
    
    # Assemble stiffness matrix.
    coordinates=inputData[5]
    meshes=inputData[3]
    numberOfMesh=len(meshes)
    materialMeshList=inputData[4]

    Iglo=[]
    Jglo=[]
    Vglo=[]

    # Integrate over non-overlapping elements.
    for i in range(len(meshes)):
        for j in range(len(meshes[i])):
            if meshes[i][j] is not None:
                coord=getCoord(coordinates,meshes[i][j])    

                if len(meshes[i][j])==3:
                    Kloc=stiffnessMatrix.triFE(coord,materialMatrix[materialMeshList[i][j]])
                elif len(meshes[i][j])==4:
                    Kloc=stiffnessMatrix.quadFE(coord,materialMatrix[materialMeshList[i][j]])
                elif len(meshes[i][j])==6:
                    Kloc=stiffnessMatrix.triQuadFE(coord,materialMatrix[materialMeshList[i][j]])
                elif len(meshes[i][j])==9:
                    Kloc=stiffnessMatrix.quadQuadFE(coord,materialMatrix[materialMeshList[i][j]])
                else: raise ValueError("Wrong element numbering!")
            
                Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,meshes[i][j])
                Iglo.extend(Iloc)
                Jglo.extend(Jloc)
                Vglo.extend(Vloc)
    
    # Integrate for the overlapping meshes.
    for i in range(len(polygons)):
        if polygons[i].triangles is None:
            # It is a traditional element.
            for j in polygons[i].incidentElements:
                if j:
                    numbering=j.numbering
                    position=j.position
                    break
            
            coord=getCoord(coordinates,numbering)
            
            if len(numbering)==3:
                Kloc=stiffnessMatrix.triFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            elif len(numbering)==4:
                Kloc=stiffnessMatrix.quadFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            elif len(numbering)==6:
                Kloc=stiffnessMatrix.triQuadFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            elif len(numbering)==9:
                Kloc=stiffnessMatrix.quadQuadFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            else: raise ValueError("Wrong element numbering!")

            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering)
            Iglo.extend(Iloc)
            Jglo.extend(Jloc)
            Vglo.extend(Vloc)
        else:
            for j in range(numberOfMesh):
                if polygons[i].incidentElements[j]:
                    numbering=polygons[i].incidentElements[j].numbering
                    position=polygons[i].incidentElements[j].position
                    coord=getCoord(coordinates,numbering)

                    Kloc=np.zeros((2*len(numbering),2*len(numbering)))

                    # Integrate the diagonal terms.
                    for k in polygons[i].triangles:
                        rho=getRho(k,position[0])
                        coordTri=coordCurvedTriangle(k,polygons[i].incidentElements,coordinates)
                        Kloc+=stiffnessMatrix.quadraticAMORE(coordTri,\
                            materialMatrix[materialMeshList[position[0]][position[1]]],coord,rho)

                    Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering)
                    Iglo.extend(Iloc)
                    Jglo.extend(Jloc)
                    Vglo.extend(Vloc)
                        
                    # Integrate the off-diagonal terms.
                    for l in range(j+1,numberOfMesh):
                        if polygons[i].incidentElements[l]:
                            numbering2=polygons[i].incidentElements[l].numbering
                            position2=polygons[i].incidentElements[l].position
                            coord2=getCoord(coordinates,numbering2)

                            Kloc=np.zeros((2*len(numbering),2*len(numbering2)))

                            for k in polygons[i].triangles:
                                rho=getRho(k,position[0])
                                rho2=getRho(k,position2[0])
                                coordTri=coordCurvedTriangle(k,polygons[i].incidentElements,coordinates)
                                Kloc+=stiffnessMatrix.quadraticAMORE(coordTri,\
                                    materialMatrix[materialMeshList[position[0]][position[1]]],coord,rho,coord2,rho2)

                            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering,numbering2)
                            Iglo.extend(Iloc)
                            Jglo.extend(Jloc)
                            Vglo.extend(Vloc)

                            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc.transpose(),numbering2,numbering)
                            Iglo.extend(Iloc)
                            Jglo.extend(Jloc)
                            Vglo.extend(Vloc)                                

    Iglo=np.array(Iglo,dtype=int)
    Jglo=np.array(Jglo,dtype=int)
    Vglo=np.array(Vglo,dtype='d')

    Kglo=scipy.sparse.coo_matrix((Vglo,(Iglo,Jglo)),shape=(2*len(coordinates),2*len(coordinates))).tocsr()
    
    print("Assembling stiffness matrix costs %s seconds."%(time.time()-startTime))
    startTime=time.time()

    # Force term.
    indForce=inputData[-2]
    if indForce[0]: # Body force is imposed.
        pass

    fglo=np.zeros((2*len(coordinates),1))

    if len(indForce)==2 or (len(indForce)==3 and indForce[2]==0):
        forceList=inputData[-1]
    elif len(indForce)==3 and indForce[2]==1:
        forceList=inputData[-1][0]
        cForceList=inputData[-1][1]

        for i in cForceList: # We assume here that no concentrated force is acting on the overlap.
            # We can remove this assumption but the code would be complicated as we need the rho functions.
            fglo[2*i[0],0]+=i[1]
            fglo[2*i[0]+1,0]+=i[2]

    nInt=3
    if indForce[1]: nInt+=2 # Customized boundary force.

    pos,wei=numericalIntegration.gaussQuad(nInt)

    nodeElements=toDC.nodeElementList(coordinates,meshes)

    for i in range(len(forceList)//2):
        node1=forceList[2*i][0]
        node2=forceList[2*i+1][0]

        force1=np.array([forceList[2*i][1:3]]).transpose()
        force2=np.array([forceList[2*i+1][1:3]]).transpose()
        force3=0.5*(force1+force2) # Linear interpolation

        # Find the element.
        elementPosition=None
        element=(set(nodeElementsPartial[node1]) & set(nodeElementsPartial[node2]))
        if element: elementPosition=element.pop()

        if elementPosition and (overlappingMesh[elementPosition[0]]\
            [elementPosition[1]].polygons[0].triangles is not None):
            # This edge is on some triangulated element.
            numbering=overlappingMesh[elementPosition[0]][elementPosition[1]].numbering
            node3=traditionalElement.findMidNode(numbering,node1,node2)

            if node3 is not None: # This edge is on some higher-order element.
                floc=np.zeros((6,1))

                if coordinates[node3]:
                    coord=np.array([coordinates[node1],coordinates[node2],coordinates[node3]])
                else:
                    coord=np.array([coordinates[node1],coordinates[node2],[0.0,0.0]])
                    coord[2,:]=0.5*(coord[0,:]+coord[1,:])

                for j in overlappingMesh[elementPosition[0]][elementPosition[1]].polygons:
                    for k in j.triangles:
                        # Test if this triangle contains the edge.
                        # If yes, integrate.
                        boundaryEdge=isTriangleContainingBoundary(k,coord[0],coord[1])
                        if boundaryEdge:
                            # In such a case, an edge of the triangle is exactly the edge defined by node1 and node2.
                            # It can be curved.
                            rho=np.array([boundaryEdge.origin.rho[elementPosition[0]],\
                                boundaryEdge.destination.rho[elementPosition[0]],0.0])
                            rho[2]=0.5*(rho[0]+rho[1])
                            coordEdge=coord

                            for l in range(nInt):
                                Nmat,Jacobian=shapeFunction.oneDQuadratic(pos[l],coordEdge)
                                rhoValue=Nmat[0,0]*rho[0]+Nmat[0,2]*rho[1]+Nmat[0,4]*rho[2]

                                if indForce[1]: # Customized boundary force
                                    xy=Nmat[0,0]*coordEdge[0]+Nmat[0,2]*coordEdge[1]+Nmat[0,4]*coordEdge[2]
                                    force=cForce.customizedForce(xy[0],xy[1])
                                else:
                                    force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3
                                
                                floc+=wei[l]*Jacobian*rhoValue*np.matmul(Nmat.transpose(),force)
                            break

                        # Test if the triangle contains part of the boundary edge.
                        boundaryEdge=isTriangleOnBoundary(k,coord[0],coord[1])
                        if boundaryEdge:
                            # It must be a straight edge.
                            rho=np.array([boundaryEdge.origin.rho[elementPosition[0]],\
                                boundaryEdge.destination.rho[elementPosition[0]]])
                            coordEdge=np.array([[boundaryEdge.origin.x,boundaryEdge.origin.y],\
                                [boundaryEdge.destination.x,boundaryEdge.destination.y]])
                            lengthEdge=lenEdge(coordEdge[0],coordEdge[1])

                            for l in range(nInt):
                                Nmat=shapeFunction.oneDLinear2(pos[l])
                                rhoValue=(Nmat@rho)[0]
                                xy=(Nmat@coordEdge).reshape(-1)
                                isoCoord=invShapeFunction.invLine(coord,xy)
                                Nmat,_=shapeFunction.oneDQuadratic(isoCoord,coord)

                                if indForce[1]: # Customized boundary force
                                    force=cForce.customizedForce(xy[0],xy[1])
                                else:
                                    force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3

                                floc+=0.5*wei[l]*lengthEdge*rhoValue*np.matmul(Nmat.transpose(),force)
                            break
                        
            else: # This edge is on some lower-order element.
                # Just copy the code for lower-order AMORE.
                floc=np.zeros((4,1))
                coord=np.array([coordinates[node1],coordinates[node2]])

                for j in overlappingMesh[elementPosition[0]][elementPosition[1]].polygons:
                    for k in j.triangles:
                        # Test if this triangle is on the boundary.
                        # If yes, integrate.
                        boundaryEdge=isTriangleOnBoundary(k,coord[0],coord[1])
                        if boundaryEdge:
                            rho=np.array([boundaryEdge.origin.rho[elementPosition[0]],\
                                boundaryEdge.destination.rho[elementPosition[0]]])
                            coordEdge=np.array([[boundaryEdge.origin.x,boundaryEdge.origin.y],\
                                [boundaryEdge.destination.x,boundaryEdge.destination.y]])
                            lengthEdge=lenEdge(coordEdge[0],coordEdge[1])
    
                            for l in range(nInt):
                                Nmat=shapeFunction.oneDLinear2(pos[l])
                                rhoValue=(Nmat@rho)[0]
                                xy=(Nmat@coordEdge).reshape(-1)
                                isoCoord=invShapeFunction.invLine(coord,xy)
                                Nmat=shapeFunction.oneDLinear(isoCoord)

                                if indForce[1]: # Customized boundary force
                                    force=cForce.customizedForce(xy[0],xy[1])
                                else:
                                    force=Nmat[0,0]*force1+Nmat[0,2]*force2

                                floc+=0.5*wei[l]*lengthEdge*rhoValue*np.matmul(Nmat.transpose(),force)
                            break

        elif elementPosition and (overlappingMesh[elementPosition[0]]\
                [elementPosition[1]].polygons[0].triangles is None):
            # A regular element.
            numbering=overlappingMesh[elementPosition[0]][elementPosition[1]].numbering
            node3=traditionalElement.findMidNode(numbering,node1,node2)
            if node3 is not None:
                floc=np.zeros((6,1))

                if coordinates[node3]:
                    coord=np.array([coordinates[node1],coordinates[node2],coordinates[node3]])
                else:
                    coord=np.array([coordinates[node1],coordinates[node2],[0.0,0.0]])
                    coord[2,:]=0.5*(coord[0,:]+coord[1,:])

                for j in range(nInt):
                    Nmat,Jacobian=shapeFunction.oneDQuadratic(pos[j],coord)

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]+Nmat[0,4]*coord[2]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3

                    floc+=wei[j]*Jacobian*np.matmul(Nmat.transpose(),force)
            else:
                floc=np.zeros((4,1))
                length=lenEdge(coordinates[node1],coordinates[node2])
                coord=np.array([coordinates[node1],coordinates[node2]])

                for j in range(nInt):
                    Nmat=shapeFunction.oneDLinear(pos[j])

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2

                    floc+=0.5*wei[j]*length*np.matmul(Nmat.transpose(),force)

        else: # This is a regular element.
            elementPosition=(set(nodeElements[node1]) & set(nodeElements[node2])).pop()
            numbering=meshes[elementPosition[0]][elementPosition[1]]
            node3=traditionalElement.findMidNode(numbering,node1,node2)
            if node3 is not None:
                floc=np.zeros((6,1))

                if coordinates[node3]:
                    coord=np.array([coordinates[node1],coordinates[node2],coordinates[node3]])
                else:
                    coord=np.array([coordinates[node1],coordinates[node2],[0.0,0.0]])
                    coord[2,:]=0.5*(coord[0,:]+coord[1,:])
                # Find the 3rd node.

                for j in range(nInt):
                    Nmat,Jacobian=shapeFunction.oneDQuadratic(pos[j],coord)

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]+Nmat[0,4]*coord[2]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3

                    floc+=wei[j]*Jacobian*np.matmul(Nmat.transpose(),force)
            else:
                floc=np.zeros((4,1))
                length=lenEdge(coordinates[node1],coordinates[node2])
                coord=np.array([coordinates[node1],coordinates[node2]])
                
                for j in range(nInt):
                    Nmat=shapeFunction.oneDLinear(pos[j])

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2

                    floc+=0.5*wei[j]*length*np.matmul(Nmat.transpose(),force)

        fglo[2*node1:2*node1+2,0]+=floc[0:2,0]
        fglo[2*node2:2*node2+2,0]+=floc[2:4,0]
        if len(floc)==6: fglo[2*node3:2*node3+2,0]+=floc[4:6,0]
    
    print("Calculating force term costs %s seconds."%(time.time()-startTime))
    startTime=time.time()
        
    # Impose constraints.
    fixList=np.zeros((2*len(coordinates),1))
    fixIndexList=np.zeros((2*len(coordinates),1),dtype=int)

    constraintList=inputData[-3]
    # Very important!!! Sort the constraints!!!
    constraintList.sort(key=lambda item:item[0])

    for i in constraintList:
        if i[1]: 
            fixList[2*i[0]]=i[3]
            fixIndexList[2*i[0]]=1
        if i[2]: 
            fixList[2*i[0]+1]=i[4]
            fixIndexList[2*i[0]+1]=1

    # Solve.
    fglo-=(Kglo.dot(fixList))

    Kglo_complete=Kglo.copy()
    Kglo=Kglo.tolil()

    count=0
    for i in constraintList:
        if i[1]:
            delete_row_lil(Kglo,2*i[0]-count)
            fglo=np.delete(fglo,2*i[0]-count)
            count+=1
        if i[2]:
            delete_row_lil(Kglo,2*i[0]+1-count)
            fglo=np.delete(fglo,2*i[0]+1-count)
            count+=1

    Kglo=Kglo.transpose()

    count=0
    for i in constraintList:
        if i[1]:
            delete_row_lil(Kglo,2*i[0]-count)
            count+=1
        if i[2]:
            delete_row_lil(Kglo,2*i[0]+1-count)
            count+=1

    print("Imposing constraints costs %s seconds."%(time.time()-startTime))
    startTime=time.time()
    
    Kglo=Kglo.tocsc()

    # factor=cholesky(Kglo)
    # disp=factor(fglo)
    disp=spsolve(Kglo,fglo)

    print("Number of equations = %s."%(len(disp)))
    print("Number of entries = %s."%(Kglo.getnnz()))
    print("Number of non-zero entries = %s."%(Kglo.count_nonzero()))

    print("Solving the linear system costs %s seconds."%(time.time()-startTime))

    # The complete displacement solution:
    displacement=np.zeros((2*len(coordinates),1))
    count=0
    for i in range(2*len(coordinates)):
        if fixIndexList[i]: 
            displacement[i]=fixList[i]
            count+=1
        else: 
            displacement[i]=disp[i-count]

    energy=0.5*displacement.transpose()@Kglo_complete@displacement

    return displacement,energy,materialMatrix

def quadraticICMAMORE(inputData,overlappingMesh,polygons,nodeElementsPartial):
    """AMORE scheme with lower-order traditional FE on the boundary meshes."""
    startTime=time.time()

    parameters=inputData[0]

    # Material matrices.
    materialList=inputData[6]

    materialMatrix=[None]*len(materialList)

    for i in range(len(materialList)):
        materialMatrix[i]=twoDMaterialMatrix(materialList[i],parameters[1])
    
    # Assemble stiffness matrix.
    coordinates=inputData[5]
    meshes=inputData[3]
    numberOfMesh=len(meshes)
    materialMeshList=inputData[4]

    Iglo=[]
    Jglo=[]
    Vglo=[]

    # Integrate over non-overlapping elements.
    for i in range(len(meshes)):
        for j in range(len(meshes[i])):
            if meshes[i][j] is not None:
                coord=getCoord(coordinates,meshes[i][j])    

                if len(meshes[i][j])==3:
                    Kloc=stiffnessMatrix.triFE(coord,materialMatrix[materialMeshList[i][j]])
                elif len(meshes[i][j])==4:
                    Kloc,_,_=stiffnessMatrix.ICMFE(coord,materialMatrix[materialMeshList[i][j]])
                elif len(meshes[i][j])==6:
                    Kloc=stiffnessMatrix.triQuadFE(coord,materialMatrix[materialMeshList[i][j]])
                elif len(meshes[i][j])==9:
                    Kloc=stiffnessMatrix.quadQuadFE(coord,materialMatrix[materialMeshList[i][j]])
                else: raise ValueError("Wrong element numbering!")
            
                Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,meshes[i][j])
                Iglo.extend(Iloc)
                Jglo.extend(Jloc)
                Vglo.extend(Vloc)
    
    # Integrate for the overlapping meshes.
    for i in range(len(polygons)):
        if polygons[i].triangles is None:
            # It is a traditional element.
            for j in polygons[i].incidentElements:
                if j:
                    numbering=j.numbering
                    position=j.position
                    break
            
            coord=getCoord(coordinates,numbering)
            
            if len(numbering)==3:
                Kloc=stiffnessMatrix.triFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            elif len(numbering)==4:
                Kloc,_,_=stiffnessMatrix.ICMFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            elif len(numbering)==6:
                Kloc=stiffnessMatrix.triQuadFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            elif len(numbering)==9:
                Kloc=stiffnessMatrix.quadQuadFE(coord,materialMatrix[materialMeshList[position[0]][position[1]]])
            else: raise ValueError("Wrong element numbering!")

            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering)
            Iglo.extend(Iloc)
            Jglo.extend(Jloc)
            Vglo.extend(Vloc)
        else:
            for j in range(numberOfMesh):
                if polygons[i].incidentElements[j]:
                    numbering=polygons[i].incidentElements[j].numbering
                    position=polygons[i].incidentElements[j].position
                    coord=getCoord(coordinates,numbering)

                    Kloc=np.zeros((2*len(numbering),2*len(numbering)))

                    # Integrate the diagonal terms.
                    for k in polygons[i].triangles:
                        rho=getRho(k,position[0])
                        coordTri=coordCurvedTriangle(k,polygons[i].incidentElements,coordinates)
                        Kloc+=stiffnessMatrix.quadraticAMORE(coordTri,\
                            materialMatrix[materialMeshList[position[0]][position[1]]],coord,rho)

                    Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering)
                    Iglo.extend(Iloc)
                    Jglo.extend(Jloc)
                    Vglo.extend(Vloc)
                        
                    # Integrate the off-diagonal terms.
                    for l in range(j+1,numberOfMesh):
                        if polygons[i].incidentElements[l]:
                            numbering2=polygons[i].incidentElements[l].numbering
                            position2=polygons[i].incidentElements[l].position
                            coord2=getCoord(coordinates,numbering2)

                            Kloc=np.zeros((2*len(numbering),2*len(numbering2)))

                            for k in polygons[i].triangles:
                                rho=getRho(k,position[0])
                                rho2=getRho(k,position2[0])
                                coordTri=coordCurvedTriangle(k,polygons[i].incidentElements,coordinates)
                                Kloc+=stiffnessMatrix.quadraticAMORE(coordTri,\
                                    materialMatrix[materialMeshList[position[0]][position[1]]],coord,rho,coord2,rho2)

                            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,numbering,numbering2)
                            Iglo.extend(Iloc)
                            Jglo.extend(Jloc)
                            Vglo.extend(Vloc)

                            Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc.transpose(),numbering2,numbering)
                            Iglo.extend(Iloc)
                            Jglo.extend(Jloc)
                            Vglo.extend(Vloc)                                

    Iglo=np.array(Iglo,dtype=int)
    Jglo=np.array(Jglo,dtype=int)
    Vglo=np.array(Vglo,dtype='d')

    Kglo=scipy.sparse.coo_matrix((Vglo,(Iglo,Jglo)),shape=(2*len(coordinates),2*len(coordinates))).tocsr()
    
    print("Assembling stiffness matrix costs %s seconds."%(time.time()-startTime))
    startTime=time.time()

    # Force term.
    indForce=inputData[-2]
    if indForce[0]: # Body force is imposed.
        pass

    fglo=np.zeros((2*len(coordinates),1))

    if len(indForce)==2 or (len(indForce)==3 and indForce[2]==0):
        forceList=inputData[-1]
    elif len(indForce)==3 and indForce[2]==1:
        forceList=inputData[-1][0]
        cForceList=inputData[-1][1]

        for i in cForceList: # We assume here that no concentrated force is acting on the overlap.
            # We can remove this assumption but the code would be complicated as we need the rho functions.
            fglo[2*i[0],0]+=i[1]
            fglo[2*i[0]+1,0]+=i[2]

    nInt=3
    if indForce[1]: nInt+=2 # Customized boundary force.

    pos,wei=numericalIntegration.gaussQuad(nInt)
    
    nodeElements=toDC.nodeElementList(coordinates,meshes)

    for i in range(len(forceList)//2):
        node1=forceList[2*i][0]
        node2=forceList[2*i+1][0]

        force1=np.array([forceList[2*i][1:3]]).transpose()
        force2=np.array([forceList[2*i+1][1:3]]).transpose()
        force3=0.5*(force1+force2) # Linear interpolation

        # Find the element.
        elementPosition=None
        element=(set(nodeElementsPartial[node1]) & set(nodeElementsPartial[node2]))
        if element: elementPosition=element.pop()

        if elementPosition and (overlappingMesh[elementPosition[0]]\
            [elementPosition[1]].polygons[0].triangles is not None):
            # This edge is on some triangulated element.
            numbering=overlappingMesh[elementPosition[0]][elementPosition[1]].numbering
            node3=traditionalElement.findMidNode(numbering,node1,node2)

            if node3 is not None: # This edge is on some higher-order element.
                floc=np.zeros((6,1))

                if coordinates[node3]:
                    coord=np.array([coordinates[node1],coordinates[node2],coordinates[node3]])
                else:
                    coord=np.array([coordinates[node1],coordinates[node2],[0.0,0.0]])
                    coord[2,:]=0.5*(coord[0,:]+coord[1,:])

                for j in overlappingMesh[elementPosition[0]][elementPosition[1]].polygons:
                    for k in j.triangles:
                        # Test if this triangle contains the edge.
                        # If yes, integrate.
                        boundaryEdge=isTriangleContainingBoundary(k,coord[0],coord[1])
                        if boundaryEdge:
                            # In such a case, an edge of the triangle is exactly the edge defined by node1 and node2.
                            # It can be curved.
                            rho=np.array([boundaryEdge.origin.rho[elementPosition[0]],\
                                boundaryEdge.destination.rho[elementPosition[0]],0.0])
                            rho[2]=0.5*(rho[0]+rho[1])
                            coordEdge=coord

                            for l in range(nInt):
                                Nmat,Jacobian=shapeFunction.oneDQuadratic(pos[l],coordEdge)
                                rhoValue=Nmat[0,0]*rho[0]+Nmat[0,2]*rho[1]+Nmat[0,4]*rho[2]

                                if indForce[1]: # Customized boundary force
                                    xy=Nmat[0,0]*coordEdge[0]+Nmat[0,2]*coordEdge[1]+Nmat[0,4]*coordEdge[2]
                                    force=cForce.customizedForce(xy[0],xy[1])
                                else:
                                    force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3
                                
                                floc+=wei[l]*Jacobian*rhoValue*np.matmul(Nmat.transpose(),force)
                            break

                        # Test if the triangle contains part of the boundary edge.
                        boundaryEdge=isTriangleOnBoundary(k,coord[0],coord[1])
                        if boundaryEdge:
                            # It must be a straight edge.
                            rho=np.array([boundaryEdge.origin.rho[elementPosition[0]],\
                                boundaryEdge.destination.rho[elementPosition[0]]])
                            coordEdge=np.array([[boundaryEdge.origin.x,boundaryEdge.origin.y],\
                                [boundaryEdge.destination.x,boundaryEdge.destination.y]])
                            lengthEdge=lenEdge(coordEdge[0],coordEdge[1])

                            for l in range(nInt):
                                Nmat=shapeFunction.oneDLinear2(pos[l])
                                rhoValue=(Nmat@rho)[0]
                                xy=(Nmat@coordEdge).reshape(-1)
                                isoCoord=invShapeFunction.invLine(coord,xy)
                                Nmat,_=shapeFunction.oneDQuadratic(isoCoord,coord)

                                if indForce[1]: # Customized boundary force
                                    force=cForce.customizedForce(xy[0],xy[1])
                                else:
                                    force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3

                                floc+=0.5*wei[l]*lengthEdge*rhoValue*np.matmul(Nmat.transpose(),force)
                            break
                        
            else: # This edge is on some lower-order element.
                # Just copy the code for lower-order AMORE.
                floc=np.zeros((4,1))
                coord=np.array([coordinates[node1],coordinates[node2]])

                for j in overlappingMesh[elementPosition[0]][elementPosition[1]].polygons:
                    for k in j.triangles:
                        # Test if this triangle is on the boundary.
                        # If yes, integrate.
                        boundaryEdge=isTriangleOnBoundary(k,coord[0],coord[1])
                        if boundaryEdge:
                            rho=np.array([boundaryEdge.origin.rho[elementPosition[0]],\
                                boundaryEdge.destination.rho[elementPosition[0]]])
                            coordEdge=np.array([[boundaryEdge.origin.x,boundaryEdge.origin.y],\
                                [boundaryEdge.destination.x,boundaryEdge.destination.y]])
                            lengthEdge=lenEdge(coordEdge[0],coordEdge[1])
    
                            for l in range(nInt):
                                Nmat=shapeFunction.oneDLinear2(pos[l])
                                rhoValue=(Nmat@rho)[0]
                                xy=(Nmat@coordEdge).reshape(-1)
                                isoCoord=invShapeFunction.invLine(coord,xy)
                                Nmat=shapeFunction.oneDLinear(isoCoord)

                                if indForce[1]: # Customized boundary force
                                    force=cForce.customizedForce(xy[0],xy[1])
                                else:
                                    force=Nmat[0,0]*force1+Nmat[0,2]*force2

                                floc+=0.5*wei[l]*lengthEdge*rhoValue*np.matmul(Nmat.transpose(),force)
                            break

        elif elementPosition and (overlappingMesh[elementPosition[0]]\
                [elementPosition[1]].polygons[0].triangles is None):
            # A regular element.
            numbering=overlappingMesh[elementPosition[0]][elementPosition[1]].numbering
            node3=traditionalElement.findMidNode(numbering,node1,node2)
            if node3 is not None:
                floc=np.zeros((6,1))

                if coordinates[node3]:
                    coord=np.array([coordinates[node1],coordinates[node2],coordinates[node3]])
                else:
                    coord=np.array([coordinates[node1],coordinates[node2],[0.0,0.0]])
                    coord[2,:]=0.5*(coord[0,:]+coord[1,:])

                for j in range(nInt):
                    Nmat,Jacobian=shapeFunction.oneDQuadratic(pos[j],coord)

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]+Nmat[0,4]*coord[2]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3

                    floc+=wei[j]*Jacobian*np.matmul(Nmat.transpose(),force)
            else:
                floc=np.zeros((4,1))
                length=lenEdge(coordinates[node1],coordinates[node2])
                coord=np.array([coordinates[node1],coordinates[node2]])

                for j in range(nInt):
                    Nmat=shapeFunction.oneDLinear(pos[j])

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2

                    floc+=0.5*wei[j]*length*np.matmul(Nmat.transpose(),force)

        else: # This is a regular element.
            elementPosition=(set(nodeElements[node1]) & set(nodeElements[node2])).pop()
            numbering=meshes[elementPosition[0]][elementPosition[1]]
            node3=traditionalElement.findMidNode(numbering,node1,node2)
            if node3 is not None:
                floc=np.zeros((6,1))

                if coordinates[node3]:
                    coord=np.array([coordinates[node1],coordinates[node2],coordinates[node3]])
                else:
                    coord=np.array([coordinates[node1],coordinates[node2],[0.0,0.0]])
                    coord[2,:]=0.5*(coord[0,:]+coord[1,:])
                # Find the 3rd node.

                for j in range(nInt):
                    Nmat,Jacobian=shapeFunction.oneDQuadratic(pos[j],coord)

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]+Nmat[0,4]*coord[2]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2+Nmat[0,4]*force3

                    floc+=wei[j]*Jacobian*np.matmul(Nmat.transpose(),force)
            else:
                floc=np.zeros((4,1))
                length=lenEdge(coordinates[node1],coordinates[node2])
                coord=np.array([coordinates[node1],coordinates[node2]])
                
                for j in range(nInt):
                    Nmat=shapeFunction.oneDLinear(pos[j])

                    if indForce[1]: # Customized boundary force
                        xy=Nmat[0,0]*coord[0]+Nmat[0,2]*coord[1]
                        force=cForce.customizedForce(xy[0],xy[1])
                    else:
                        force=Nmat[0,0]*force1+Nmat[0,2]*force2

                    floc+=0.5*wei[j]*length*np.matmul(Nmat.transpose(),force)

        fglo[2*node1:2*node1+2,0]+=floc[0:2,0]
        fglo[2*node2:2*node2+2,0]+=floc[2:4,0]
        if len(floc)==6: fglo[2*node3:2*node3+2,0]+=floc[4:6,0]
    
    print("Calculating force term costs %s seconds."%(time.time()-startTime))
    startTime=time.time()
        
    # Impose constraints.
    fixList=np.zeros((2*len(coordinates),1))
    fixIndexList=np.zeros((2*len(coordinates),1),dtype=int)

    constraintList=inputData[-3]
    # Very important!!! Sort the constraints!!!
    constraintList.sort(key=lambda item:item[0])

    for i in constraintList:
        if i[1]: 
            fixList[2*i[0]]=i[3]
            fixIndexList[2*i[0]]=1
        if i[2]: 
            fixList[2*i[0]+1]=i[4]
            fixIndexList[2*i[0]+1]=1

    # Solve.
    fglo-=(Kglo.dot(fixList))

    Kglo_complete=Kglo.copy()
    Kglo=Kglo.tolil()

    count=0
    for i in constraintList:
        if i[1]:
            delete_row_lil(Kglo,2*i[0]-count)
            fglo=np.delete(fglo,2*i[0]-count)
            count+=1
        if i[2]:
            delete_row_lil(Kglo,2*i[0]+1-count)
            fglo=np.delete(fglo,2*i[0]+1-count)
            count+=1

    Kglo=Kglo.transpose()

    count=0
    for i in constraintList:
        if i[1]:
            delete_row_lil(Kglo,2*i[0]-count)
            count+=1
        if i[2]:
            delete_row_lil(Kglo,2*i[0]+1-count)
            count+=1

    print("Imposing constraints costs %s seconds."%(time.time()-startTime))
    startTime=time.time()
    
    Kglo=Kglo.tocsc()

    # factor=cholesky(Kglo)
    # disp=factor(fglo)
    disp=spsolve(Kglo,fglo)

    print("Number of equations = %s."%(len(disp)))
    print("Number of stored entries = %s."%(Kglo.getnnz()))
    print("Number of non-zero entries = %s."%(Kglo.count_nonzero()))

    print("Solving the linear system costs %s seconds."%(time.time()-startTime))

    # The complete displacement solution:
    displacement=np.zeros((2*len(coordinates),1))
    count=0
    for i in range(2*len(coordinates)):
        if fixIndexList[i]: 
            displacement[i]=fixList[i]
            count+=1
        else: 
            displacement[i]=disp[i-count]

    energy=0.5*displacement.transpose()@Kglo_complete@displacement

    return displacement,energy,materialMatrix

def twoDMaterialMatrix(material,problemType):
    """Input:
    
    material: [E,nu];
    
    problemType: 1 -- Plane stress; 2 -- Plane strain."""

    if problemType==1: # Plane stress
        Dmat=np.array([[material[0]/(1.0-material[1]**2),material[0]*material[1]/(1.0-material[1]**2),0.0],\
                       [material[0]*material[1]/(1.0-material[1]**2),material[0]/(1.0-material[1]**2),0.0],\
                       [0.0,0.0,material[0]/2.0/(1.0+material[1])]])
    elif problemType==2: #Plane strain
        cc=material[0]*(1.0-material[1])/(1.0+material[1])/(1.0-2.0*material[1])
        Dmat=np.array([[cc,cc*material[1]/(1.0-material[1]),0.0],\
                       [cc*material[1]/(1.0-material[1]),cc,0.0],\
                       [0.0,0.0,cc*(1.0-2.0*material[1])/(2.0*(1.0-material[1]))]])
    else: raise ValueError("No such a problem type!")
    
    return Dmat

def coordCurvedTriangle(HalfEdge,incidentElements,coordinates):
    edge=HalfEdge
    startingEdge=edge

    for i in incidentElements:
        if i and len(i.numbering)>4: # A higher-order element
            edge=startingEdge

            if len(i.numbering)==9: 
                numberOfNode=4
                temp=[1,2,3,0]
            elif len(i.numbering)==6: 
                numberOfNode=3
                temp=[1,2,0]
            else: raise ValueError

            # Find the curved edge in i.
            curvedEdge=[]
            for j in range(numberOfNode):
                if coordinates[i.numbering[j+numberOfNode]] is not None: curvedEdge.append(j)
            
            ind=[0,0,0]
            ind2=[None,None,None]
            for j in range(3):
                for k in curvedEdge:
                    if isEqual(edge.origin,coordinates[i.numbering[k]]) and \
                        isEqual(edge.destination,coordinates[i.numbering[temp[k]]]):
                        ind[j]=1
                        ind2[j]=k
                        break
                edge=edge.next

            if sum(ind)>0:
                edge=startingEdge

                coord=np.zeros((6,2))

                for j in range(3):
                    coord[j,:]=np.array([edge.origin.x,edge.origin.y])

                    if ind2[j] is not None:
                        coord[j+3,:]=np.array(coordinates[i.numbering[ind2[j]+numberOfNode]])
                    else:
                        coord[j+3,:]=np.array(midEdge(edge))

                    edge=edge.next
                
                return coord

        elif i and len(i.numbering)<=4: return stiffnessMatrix.getCoordFromHalfEdge(startingEdge)

    return stiffnessMatrix.getCoordFromHalfEdge(startingEdge)

def midEdge(edge):
    return [(edge.origin.x+edge.destination.x)/2,(edge.origin.y+edge.destination.y)/2]

def isEqual(vertex,coordinate):
    if vertex.x==coordinate[0] and vertex.y==coordinate[1]: return True
    else: return False

def getCoord(coordinates,numbering):
    coord=np.zeros((len(numbering),2))

    if len(numbering)==3 or len(numbering)==4: # 1st-order element.
        for i in range(len(numbering)):
            coord[i,:]=np.array(coordinates[numbering[i]])
    elif len(numbering)==6: # 2nd-order triangle.
        for i in range(3):
            coord[i,:]=np.array(coordinates[numbering[i]])
        temp=[1,2,0]
        for i in range(3,6):
            if coordinates[numbering[i]]: coord[i,:]=np.array(coordinates[numbering[i]])
            else:
                coord[i,:]=0.5*(np.array(coordinates[numbering[i-3]])+np.array(coordinates[numbering[temp[i-3]]]))
    elif len(numbering)==9: # 2nd-order quad.
        for i in range(4):
            coord[i,:]=np.array(coordinates[numbering[i]])
        temp=[1,2,3,0]
        for i in range(4,8):
            if coordinates[numbering[i]]: coord[i,:]=np.array(coordinates[numbering[i]])
            else:
                coord[i,:]=0.5*(np.array(coordinates[numbering[i-4]])+np.array(coordinates[numbering[temp[i-4]]]))
        if coordinates[numbering[8]]: coord[8,:]=np.array(coordinates[numbering[8]])
        else:
            coord[8,:]=0.25*(np.array(coordinates[numbering[0]])+\
                            np.array(coordinates[numbering[1]])+\
                            np.array(coordinates[numbering[2]])+\
                            np.array(coordinates[numbering[3]]))
    else: raise ValueError

    return coord

def lenEdge(coord1,coord2):
    return ((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2)**0.5

def delete_row_lil(matrix,i):
    if not isinstance(matrix,scipy.sparse.lil_matrix):
        raise ValueError("The matrix should be in LIL format!")
    matrix.rows=np.delete(matrix.rows,i)
    matrix.data=np.delete(matrix.data,i)
    matrix._shape=(matrix._shape[0]-1,matrix._shape[1])

def getRho(edge,i): # This is actually the weight functions.
    rho=np.array([edge.origin.rho[i],edge.next.origin.rho[i],edge.next.next.origin.rho[i]])
    return rho

def isTriangleOnBoundary(edge,coord1,coord2):
    numberOfEgde=triangulation.countEdges(edge)

    length=math.sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2)
    length1=math.sqrt((edge.origin.x-coord2[0])**2+(edge.origin.y-coord2[1])**2)

    for _ in range(numberOfEgde):
        length2=math.sqrt((edge.destination.x-coord2[0])**2+(edge.destination.y-coord2[1])**2)
        
        if abs(planeSweepMeshes.crossProductBasic(edge.origin.x-coord2[0],\
            edge.origin.y-coord2[1],coord1[0]-coord2[0],coord1[1]-coord2[1]))<=10**(-13)*length*length1\
                and \
            abs(planeSweepMeshes.crossProductBasic(edge.destination.x-coord2[0],\
            edge.destination.y-coord2[1],coord1[0]-coord2[0],coord1[1]-coord2[1]))<=10**(-13)*length*length2:
            return edge
        length1=length2
        edge=edge.next

    return None

def isTriangleContainingBoundary(edge,coord1,coord2):
    numberOfEgde=triangulation.countEdges(edge)

    for _ in range(numberOfEgde):
        if (edge.origin.x==coord1[0] and edge.origin.y==coord1[1] and \
            edge.destination.x==coord2[0] and edge.destination.y==coord2[1]):
            return edge
        elif (edge.origin.x==coord2[0] and edge.origin.y==coord2[1] and \
                    edge.destination.x==coord1[0] and edge.destination.y==coord1[1]):
            return edge.twin
    
        edge=edge.next

    return None
