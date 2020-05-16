import sys

def twoDListMin(elements):
    """Return the min of a 2D list.

    Input: A list of lists of numbers. It can not be empty."""
    
    if len(elements)>1:
        theMin=min(twoDListMin(elements[0 : (len(elements)//2)]),\
            twoDListMin(elements[(len(elements)//2) : len(elements)]))
    else:
        if len(elements)==0: return 0
        theMin=min(elements[0])

    return theMin

def strNoComment(string,comment=None):
    """Remove the comments from a line string:
    If the string starts with the comment symbol, return None.
    Otherwise, the chars after the first comment symbol are removed."""
    string=string.strip() # Remove the white space on both sides.
    if comment is not None: 
        sList=string.split(comment)
        return sList[0]
    else: return string

def skipSpace(file,comment):
    """Fast move to the next non-empty line.
    Return the line as a string without any comment."""
    count=0 # Count how many rows are skipped.

    a=strNoComment(file.readline(),comment)
    while len(a)==0: # The line just read is empty. Read the next line.
        a=strNoComment(file.readline(),comment)
        count+=1
        assert (count<=1000),"Invalid input leads to looping!"
    return a

def txtInput(path):
    """Given a .txt file name in the inputFiles folder, read data from it
    for FE analysis."""
    # Open the file as readable
    try:
        file=open("./inputFiles/"+path+".txt",'r')
    except FileNotFoundError:
        raise
    except:
        print("Unexpected Error")
        sys.exit()

    # We skip all the problem description.
    a=skipSpace(file,"#").replace(',',' ').split()
    # Now we find the first line that indicates the mesh parameters.
    # We can also use "," to split different numbers.
    
    assert (len(a)==2),"Error in Number of Parameters"
    parameters=[int(i) for i in a]

    a=skipSpace(file,"#").replace(',',' ').split()

    # If numSolver==6|7 (OFE|AMORE), we need more information.
    basis=indFEwC=None
    if parameters[0]==6 or parameters[0]==7:
        assert (len(a)==2),"Error in Number of Parameters"
        basis=int(a[0])
        indFEwC=int(a[1])
        a=skipSpace(file,"#").replace(',',' ').split()

    # We now construct the element connectivity matrix and the element matrial list.
    numberOfMesh=int(a[0])
    del a[0]
    numberOfElement=[int(j) for j in a]

    a=skipSpace(file,"#").replace(',',' ').split()

    meshes=[None]*numberOfMesh
    materialMeshList=[None]*numberOfMesh

    for j in range(numberOfMesh):
        elements=[]
        eleMaterial=[]
        
        for i in range(numberOfElement[j]):
            element=[int(k) for k in a] # The element numbering should be counterclockwise.
            elements.append(element[0:len(element)-1])
            eleMaterial.append(element[-1])
            a=skipSpace(file,"#").replace(',',' ').split()
        
        meshes[j]=elements
        materialMeshList[j]=eleMaterial
    
    totalMesh=[None]*numberOfMesh
    for i in range(numberOfMesh):
        totalMesh[i]=twoDListMin(meshes[i])

    assert(min(totalMesh)==0),"The node numbering should start from 0."
        
    # Now we collect the coordinates.
    numberOfNode=int(a[0])
    a=skipSpace(file,"#").replace(',',' ').split()

    coordinates=[]

    for i in range(numberOfNode):
        try:
            coordinate=[float(j) for j in a]
        except:
            coordinate=None # In such a case, this is a mid-edge node.
        coordinates.append(coordinate)
        a=skipSpace(file,"#").replace(',',' ').split()

    numberOfMaterial=int(a[0])
    a=skipSpace(file,"#").replace(',',' ').split()

    materialList=[]
    # Material properties
    for i in range(numberOfMaterial):
        material=[float(i) for i in a]
        materialList.append(material)
        a=skipSpace(file,"#").replace(',',' ').split()

    numberOfConstraint=int(a[0])
    a=skipSpace(file,"#").replace(',',' ').split()

    constraintList=[]
    # Constraints
    for i in range(numberOfConstraint):
        if len(a)>5:
            constraint=[int(a[0]),int(a[1]),int(a[2]),float(a[3]),float(a[4]),float(a[5]),\
                        float(a[6]),int(a[7]),int(a[8])]
        else:
            constraint=[int(a[0]),int(a[1]),int(a[2]),float(a[3]),float(a[4])]
        constraintList.append(constraint)
        a=skipSpace(file,"#").replace(',',' ').split()

    if parameters[0]!=6 and parameters[0]!=7: # Traditional elements
        for i in range(numberOfConstraint):
            constraintList[i]=constraintList[i][0:5]

    # Force terms
    indForce=[int(i) for i in a] # Body force & Customized boundary force & [Concentrated force]

    if len(indForce)==2 or (len(indForce)==3 and indForce[2]==0): # No concentrated force
        # Number of force pairs
        a=skipSpace(file,"#").replace(',',' ').split()
        assert (len(a)==1),"Invalid input on boundary forces!"

        numberOfForcePair=int(a[0])

        if numberOfForcePair!=0:
            a=skipSpace(file,"#").replace(',',' ').split()

            assert (len(a)==3),"Error in boundary tractions!"
        forcePairList=[]
        
        for i in range(numberOfForcePair):
            forcePair=[int(a[0]),float(a[1]),float(a[2])]
            forcePairList.append(forcePair)
            a=skipSpace(file,"#").replace(',',' ').split()
            forcePair=[int(a[0]),float(a[1]),float(a[2])]
            forcePairList.append(forcePair)
            if i!=numberOfForcePair-1:
                a=skipSpace(file,"#").replace(',',' ').split()

        return (parameters,basis,indFEwC,meshes,materialMeshList,coordinates,materialList,\
                constraintList,indForce,forcePairList)

    elif len(indForce)==3 and indForce[2]==1: # In such a case, indForce[2] is indConcentrated.
        # Number of concentrated forces
        a=skipSpace(file,"#").replace(',',' ').split() # In this case, len(a) = 2.
        assert (len(a)==2),"Invalid input on boundary forces!"

        numberOfForcePair=int(a[0])
        numberOfCForce=int(a[1])

        if numberOfForcePair!=0:
            a=skipSpace(file,"#").replace(',',' ').split()

            assert (len(a)==3),"Error in boundary tractions!"
        forcePairList=[]
        
        for i in range(numberOfForcePair):
            forcePair=[int(a[0]),float(a[1]),float(a[2])]
            forcePairList.append(forcePair)
            a=skipSpace(file,"#").replace(',',' ').split()
            forcePair=[int(a[0]),float(a[1]),float(a[2])]
            forcePairList.append(forcePair)
            if i!=numberOfForcePair-1:
                a=skipSpace(file,"#").replace(',',' ').split()

        if numberOfCForce!=0:
            a=skipSpace(file,"#").replace(',',' ').split()

            assert (len(a)==3),"Error in boundary tractions!"
        cForceList=[]
        
        for i in range(numberOfCForce):
            cForce=[int(a[0]),float(a[1]),float(a[2])]
            cForceList.append(cForce)

            if i!=numberOfCForce-1:
                a=skipSpace(file,"#").replace(',',' ').split()

        return (parameters,basis,indFEwC,meshes,materialMeshList,coordinates,materialList,\
                constraintList,indForce,[forcePairList,cForceList]) # Now the output is different.

    else: raise ValueError

if __name__=="__main__":
    pass