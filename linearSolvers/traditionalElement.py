import numpy as np 
import scipy
# from sksparse.cholmod import cholesky # It works (from Terminal).
from scipy.sparse.linalg import spsolve
import time
import sys, os
sys.path.append(os.path.dirname(sys.path[0]))
from elementLibrary import stiffnessMatrix, shapeFunction
from otherFunctions import numericalIntegration
from linearSolvers import AMORE
from meshTools import toDC

def lowerOrderFE(inputData):
    """The solver for lower order finite elements (4-node quads & 3-node triangles)."""
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
    materialMeshList=inputData[4]

    Iglo=[]
    Jglo=[]
    Vglo=[]

    for i in range(len(meshes[0])):
        coord=AMORE.getCoord(coordinates,meshes[0][i])

        if len(meshes[0][i])==3:
            Kloc=stiffnessMatrix.triFE(coord,materialMatrix[materialMeshList[0][i]])
        elif len(meshes[0][i])==4:
            Kloc=stiffnessMatrix.quadFE(coord,materialMatrix[materialMeshList[0][i]])
        else: raise ValueError("Wrong element numbering!")
        
        Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,meshes[0][i])
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

    forceList=inputData[-1]

    nInt=2
    if indForce[1]: nInt+=3 # Customized boundary force.

    pos,wei=numericalIntegration.gaussQuad(nInt)
    fglo=np.zeros((2*len(coordinates),1))

    for i in range(len(forceList)//2):
        node1=forceList[2*i][0]
        node2=forceList[2*i+1][0]
        length=lenEdge(coordinates[node1],coordinates[node2])

        force1=np.array([forceList[2*i][1:3]]).transpose()
        force2=np.array([forceList[2*i+1][1:3]]).transpose()

        floc=np.zeros((4,1))

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

def ICMFE(inputData):
    """The solver for (4-node) ICM finite elements. Warning: The code is only for squares. 
    For general quadrilaterals, the formulation needs to be modified to pass patch tests."""
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
    materialMeshList=inputData[4]

    Iglo=[]
    Jglo=[]
    Vglo=[]

    for i in range(len(meshes[0])):
        coord=AMORE.getCoord(coordinates,meshes[0][i])

        if len(meshes[0][i])==4:
            Kloc,_,_=stiffnessMatrix.ICMFE(coord,materialMatrix[materialMeshList[0][i]])
        else: raise ValueError("Wrong element numbering!")
        
        Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,meshes[0][i])
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

    forceList=inputData[-1]

    nInt=2
    if indForce[1]: nInt+=3 # Customized boundary force.

    pos,wei=numericalIntegration.gaussQuad(nInt)
    fglo=np.zeros((2*len(coordinates),1))

    for i in range(len(forceList)//2):
        node1=forceList[2*i][0]
        node2=forceList[2*i+1][0]
        length=lenEdge(coordinates[node1],coordinates[node2])

        force1=np.array([forceList[2*i][1:3]]).transpose()
        force2=np.array([forceList[2*i+1][1:3]]).transpose()

        floc=np.zeros((4,1))

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
    print("Number of non-zero sparse matrix entries = %s."%Kglo.count_nonzero())

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

def quadraticFE(inputData):
    """The solver for second order finite elements (9-node quads & 6-node triangles)."""
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
    materialMeshList=inputData[4]

    nodeElements=toDC.nodeElementList(coordinates,meshes)

    Iglo=[]
    Jglo=[]
    Vglo=[]

    for i in range(len(meshes[0])):
        coord=AMORE.getCoord(coordinates,meshes[0][i])

        if len(meshes[0][i])==6:
            Kloc=stiffnessMatrix.triQuadFE(coord,materialMatrix[materialMeshList[0][i]])
        elif len(meshes[0][i])==9:
            Kloc=stiffnessMatrix.quadQuadFE(coord,materialMatrix[materialMeshList[0][i]])
        else: raise ValueError("Wrong element numbering!")
        
        Iloc,Jloc,Vloc=stiffnessMatrix.sparsifyElementMatrix(Kloc,meshes[0][i])
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

    forceList=inputData[-1]

    nInt=3
    if indForce[1]: nInt+=2 # Customized boundary force.

    pos,wei=numericalIntegration.gaussQuad(nInt)
    fglo=np.zeros((2*len(coordinates),1))

    for i in range(len(forceList)//2):
        node1=forceList[2*i][0]
        node2=forceList[2*i+1][0]
        # length=lenEdge(coordinates[node1],coordinates[node2])

        # Find the element.
        elementPosition=(set(nodeElements[node1]) & set(nodeElements[node2])).pop()
        numbering=meshes[elementPosition[0]][elementPosition[1]]
        node3=findMidNode(numbering,node1,node2)

        force1=np.array([forceList[2*i][1:3]]).transpose()
        force2=np.array([forceList[2*i+1][1:3]]).transpose()

        floc=np.zeros((6,1))
        coord=np.array([coordinates[node1],coordinates[node2],[0.0,0.0]])
        if coordinates[node3]: coord[2,:]=np.array(coordinates[node3])
        else: coord[2,:]=0.5*(coord[0,:]+coord[1,:])

        for j in range(nInt):
            # Only support linear force distribution.
            # Otherwise, use customized boundary force.
            Nmat=shapeFunction.oneDLinear(pos[j])
            force=Nmat[0,0]*force1+Nmat[0,2]*force2

            quadNmat,Jacobian=shapeFunction.oneDQuadratic(pos[j],coord)

            floc+=wei[j]*Jacobian*np.matmul(quadNmat.transpose(),force)

        fglo[2*node1:2*node1+2,0]+=floc[0:2,0]
        fglo[2*node2:2*node2+2,0]+=floc[2:4,0]
        fglo[2*node3:2*node3+2,0]+=floc[4:6,0]
    
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

def findMidNode(numbering,node1,node2):
    if len(numbering)==6:
        temp=[1,2,0]
        for i in range(3):
            if (numbering[i]==node1 and numbering[temp[i]]==node2) or \
                (numbering[i]==node2 and numbering[temp[i]]==node1):
                return numbering[i+3]
    elif len(numbering)==9:
        temp=[1,2,3,0]
        for i in range(4):
            if (numbering[i]==node1 and numbering[temp[i]]==node2) or \
                (numbering[i]==node2 and numbering[temp[i]]==node1):
                return numbering[i+4]

    return None

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
    else: raise ValueError("No such problem type!")
    
    return Dmat

def lenEdge(coord1,coord2):
    return ((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2)**0.5

def delete_row_lil(matrix,i):
    if not isinstance(matrix,scipy.sparse.lil_matrix):
        raise ValueError("The matrix should be in LIL format!")
    matrix.rows=np.delete(matrix.rows,i)
    matrix.data=np.delete(matrix.data,i)
    matrix._shape=(matrix._shape[0]-1,matrix._shape[1])

if __name__=="__main__":
    pass