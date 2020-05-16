import sys,os
sys.path.append(os.path.dirname(sys.path[0]))
from computGeometry import planeSweepMeshes

def nodeElementList(coordinates,meshes): # Total connectivity list
    """Construct the node-element inverse connectivity"""
    nodeElements=[None]*len(coordinates)
    for i in range(len(coordinates)):
        nodeElements[i]=[]

    for i in range(len(meshes)):
        for j in range(len(meshes[i])):
            if meshes[i][j] is not None:
                if len(meshes[i][j])==9 or len(meshes[i][j])==4:
                    numberOfNode=4
                elif len(meshes[i][j])==6 or len(meshes[i][j])==3:
                    numberOfNode=3
                else: raise ValueError

                for k in range(numberOfNode):
                    nodeElements[meshes[i][j][k]].append((i,j)) 
                    # List is not hashable (for set operations), so we use tuple here.

    return nodeElements

def nodeElementPartialList(coordinates,overlappingMesh): 
    # Partial connectivity list for overlapping elements
    nodeElements=[None]*len(coordinates)
    for i in range(len(coordinates)):
        nodeElements[i]=[]

    for i in range(len(overlappingMesh)):
        for j in range(len(overlappingMesh[i])):
            if len(overlappingMesh[i][j].numbering)==9 or len(overlappingMesh[i][j].numbering)==4:
                numberOfNode=4
            elif len(overlappingMesh[i][j].numbering)==6 or len(overlappingMesh[i][j].numbering)==3:
                numberOfNode=3
            else: raise ValueError

            for k in range(numberOfNode):
                # We only need to deal with these main nodes.
                nodeElements[overlappingMesh[i][j].numbering[k]].append((i,j))

    return nodeElements

######

class face(object):
    def __init__(self,element,outer):
        self.element=element
        self.outer=outer
        self.label=None # Indicating which mesh the element belongs to
        self.incidentElements=None
        self.triangles=None

    def __str__(self):
        return str(self.element)

class halfEdge(object):
    def __init__(self,origin,face):
        self.origin=origin
        self.destination=None
        self.incidentFace=face
        self.twin=None
        self.next=None
        self.prev=None
        self.line=None
        self.faceList=[] # Store the faces that contain the left side of this half edge.
        self.newFace=None # After establishing the connections between half-edges, we create a new face for each subregion.

    def __str__(self):
        return str(self.line)
######

######
# Find the next edge
def findNext(edge,edgeList):
    """Find the next edge from the given list."""
    for i in edgeList:
        if i.prev is None and edge.destination is i.origin: return i
    return None

# Find the previous edge
def findPrev(edge,edgeList):
    """Find the previous edge from the given list."""
    for i in edgeList:
        if i.next is None and edge.origin is i.destination: return i
    return None
######

def toDoublyConnected(coordinates,overlappingMesh):
    """Convert a FE mesh to its doubly-connected edge list representation"""
    ######
    # coordinates: 2D list
    # elements: 2D list
    ######
    # We require now that the node numbering starts from 0.

    ######
    # Initiate lists
    vertexList=[]
    halfEdgeList=[]
    faceList=[]
    lineList=[]
    ######

    ######
    
    nodeElementsPartial=nodeElementPartialList(coordinates,overlappingMesh)
        
    # Convert vertex list
    for i in range(len(coordinates)):
        if coordinates[i] is not None:
            vertexList.append(planeSweepMeshes.point_2D(coordinates[i][0],coordinates[i][1],i))
            vertexList[-1].overlappingVertex=[None]*len(overlappingMesh) 
            # In the worst case, each vertex may have len(overlappingMesh) overlapping vertices, including itself.
    
            if len(nodeElementsPartial[i])>0: 
                vertexList[i].label=nodeElementsPartial[i][0][0]
                # If this node is in overlapping elements, give it the label indicating which mesh it lies in.
                # Elements in different meshes can not share a common node, though the coordinates may be coincident.
                vertexList[-1].overlappingVertex[vertexList[-1].label]=vertexList[-1]
                # The overlapping vertex in the labelled mesh is of course the vertex itself.

        else: vertexList.append(None)

    # Half edge list and face list
    for i in range(len(overlappingMesh)):
        for j in range(len(overlappingMesh[i])):
            if len(overlappingMesh[i][j].numbering)==9 or len(overlappingMesh[i][j].numbering)==4:
                numberOfNode=4
            elif len(overlappingMesh[i][j].numbering)==6 or len(overlappingMesh[i][j].numbering)==3:
                numberOfNode=3
            else: raise ValueError

            for k in overlappingMesh[i][j].numbering[0:numberOfNode]:
                overlappingMesh[i][j].vertex.append(vertexList[k])

            faceList.append(face(overlappingMesh[i][j],None))
            faceList[-1].label=i
            overlappingMesh[i][j].face=faceList[-1]
            
            if numberOfNode==3: 
                temp=[1,2,0]
                temp1=[2,0,1]
            else: 
                temp=[1,2,3,0]
                temp1=[3,0,1,2]

            for k in range(numberOfNode):
                halfEdgeList.append(halfEdge(vertexList[overlappingMesh[i][j].numbering[k]],faceList[-1]))
                halfEdgeList[-1].destination=vertexList[overlappingMesh[i][j].numbering[temp[k]]]

                if vertexList[overlappingMesh[i][j].numbering[k]].incidentEdge is None:
                    vertexList[overlappingMesh[i][j].numbering[k]].incidentEdge=halfEdgeList[-1]
            
            faceList[-1].outer=halfEdgeList[-1]

            for k in range(numberOfNode):
                halfEdgeList[k-numberOfNode].next=halfEdgeList[temp[k]-numberOfNode]
                halfEdgeList[k-numberOfNode].prev=halfEdgeList[temp1[k]-numberOfNode]

    # Now we assign twin for each edge.
    for i in halfEdgeList:
        if i.twin is None:
            temp=list(set(nodeElementsPartial[i.origin.num]) & set(nodeElementsPartial[i.destination.num]))
            if len(temp)==2: # At most two elements can share the same edge.
                if i.incidentFace.element is overlappingMesh[temp[0][0]][temp[0][1]]:
                    otherElement=temp[1]
                else:
                    otherElement=temp[0]

                tempEdge=overlappingMesh[otherElement[0]][otherElement[1]].face.outer
                # Found an edge in the other element.

                for j in range(len(overlappingMesh[otherElement[0]][otherElement[1]].numbering)):
                    # Loop to find the exact twin.
                    if tempEdge.origin is i.destination:
                        i.twin=tempEdge
                        tempEdge.twin=i
                        break
                    else:
                        tempEdge=tempEdge.next

    # The outer boundary
    lenEdge=len(halfEdgeList) # Number of regular edges (established directly from the element numberings)

    for i in range(lenEdge):
        if halfEdgeList[i].twin is None: # This must be on the outer loop.
            halfEdgeList.append(halfEdge(halfEdgeList[i].destination,None)) # The incident face is None.
            halfEdgeList[-1].destination=halfEdgeList[i].origin
            halfEdgeList[-1].twin=halfEdgeList[i]
            halfEdgeList[i].twin=halfEdgeList[-1]

    for i in range(lenEdge,len(halfEdgeList)):
        if halfEdgeList[i].next is None:
            nextEdge=findNext(halfEdgeList[i],halfEdgeList[i+1:len(halfEdgeList)])
            halfEdgeList[i].next=nextEdge
            nextEdge.prev=halfEdgeList[i]

        if halfEdgeList[i].prev is None:
            prevEdge=findPrev(halfEdgeList[i],halfEdgeList[i+1:len(halfEdgeList)])
            halfEdgeList[i].prev=prevEdge
            prevEdge.next=halfEdgeList[i]
    ######

    # We calculate the line list.
    for i in halfEdgeList:
        if i.line is None:
            lineList.append(planeSweepMeshes.line(i.origin,i.destination))
            i.line=lineList[-1]
            i.twin.line=lineList[-1]
            lineList[-1].halfEdge=[i,i.twin]
            i.line.upperPoint.incidentU.append(i.line)
            i.line.lowerPoint.incidentL.append(i.line)
            # incidentL and incidentU have now been correctly established for each mesh.

    ######
    # Remove unused vertices.
    newVertexList=[]
    for i in range(len(coordinates)):
        if nodeElementsPartial[i]: 
            vertexList[i].weight=1.0
            newVertexList.append(vertexList[i])

    del vertexList

    return (newVertexList,faceList,halfEdgeList,lineList,nodeElementsPartial)

def printPolygon(polygon,coordinates):
    startingEdge=polygon.outer
    edge=startingEdge

    coord=[]
    coord.append(str(edge.origin))
    edge=edge.next
    while edge is not startingEdge:
        coord.append(str(edge.origin))
        edge=edge.next

    print('Coordinates of this polygon is %s.'%coord)

    for i in range(len(polygon.incidentElements)):
        coord=[]
        if polygon.incidentElements[i]:
            if len(polygon.incidentElements[i].numbering)==4 or len(polygon.incidentElements[i].numbering)==9:
                numberOfNode=4
            else: numberOfNode=3
            for j in range(numberOfNode):
                coord.append(coordinates[polygon.incidentElements[i].numbering[j]])
            print('The No. %s incident element has coordinates %s.'%(i,coord))
        else: print('The No. %s incident element is None.'%i)

def printTriangle(halfEdge):
    print("Coordinates of the triangle: ")
    edge=halfEdge
    for _ in range(3):
        print(edge.origin)
        print('rho = %s'%edge.origin.rho)
        edge=edge.next

######
if __name__=="__main__":
    print()