# The plane sweep algorithm for overlapping meshes.
import math
import sys, os
import numpy as np 
sys.path.append(os.path.dirname(sys.path[0]))
from computGeometry import triangulation
from DataStructure import bstAVL,bstAVLLines
from inputFunctions import inputF
from plotTools import pM
from meshTools import toDC
from elementLibrary import invShapeFunction

global scale, tol 
scale=[1,1]
tol=10**(-14)

class point_2D: # The vertex class.
    def __init__(self,x_coord,y_coord,number=None):
        self.num=number
        self.overlappingVertex=None
        self.label=None
        self.weight=None # We use it to establish partition of unity functions.
        self.rho=None # Value of the weight functions.
        self.x=np.longdouble(x_coord)
        self.y=np.longdouble(y_coord)
        self.incidentU=[] # All lines that have self as the upper end point.
        self.incidentL=[] # All lines that have self as the lower end point.
        self.incidentC=[] # All lines that have self as an interior point.
        self.leftList=None # We need to store the immediate left edge in each mesh to the event point.
                           # This will be used to determine what elements the new face lies in.
                           # If the mesh is involved, we store one incident half-edge originating from the event point.
        self.incidentEdge=None
        self.involvedMesh=None # We maintain a list that gives which meshes are involved at the event point.
    
    def __gt__(self,other):
        if self.y != other.y:
            return self.y > other.y
        else:
            return self.x < other.x # This special order is given in the book of de Berg et al.
    
    def __lt__(self,other):
        if self.y != other.y:
            return self.y < other.y
        else:
            return self.x > other.x
    
    def __eq__(self,other):
        if abs(self.x-other.x)<tol*scale[0] and abs(self.y-other.y)<tol*scale[1]:
            return True
        else:
            return False
    
    def __ge__(self,other):
        return self>other or self==other
    
    def __le__(self,other):
        return self<other or self==other

    def __ne__(self,other):
        return not self==other

    def __str__(self):
        return str([self.x,self.y])

class line:
    def __init__(self,point1,point2):
        if point1>point2: # Put the upper point on top.
            self.upperPoint=point1
            self.lowerPoint=point2
        else:
            self.upperPoint=point2
            self.lowerPoint=point1

        self.halfEdge=None # Two half edges corresponding to this line.
    
    def length(self):
        return math.sqrt((self.upperPoint.x-self.lowerPoint.x)**2+\
            (self.upperPoint.y-self.lowerPoint.y)**2)

    def __str__(self):
        return str(self.upperPoint)+' -- '+str(self.lowerPoint)

def isOnLine(line,point):
    """Test if the point is on the line segment."""
    deltaX=line.upperPoint.x-line.lowerPoint.x 
    deltaY=line.lowerPoint.y-line.upperPoint.y

    tol1 = (abs(deltaY) * scale[0] + abs(deltaX) * scale[1]) * tol # Tolerance

    if abs((point.x-line.lowerPoint.x)*deltaY+(point.y-line.lowerPoint.y)*deltaX)<=tol1 and \
        (line.lowerPoint.y-scale[1]*tol<=point.y<=line.upperPoint.y+scale[1]*tol) and \
            (point.x-line.lowerPoint.x)*(point.x-line.upperPoint.x)<=scale[0]*tol*(abs(deltaX)+scale[0]*tol):
        return True
    else: return False

def calHorizonInter(line,point):
    """Calculate the intersection of a given line segment and the horizontal line through the given point."""

    if isOnLine(line,point): return point
    elif line.upperPoint.y==line.lowerPoint.y:
        # In this case, the edge has not yet been added into the status or has already left.
        # We can return every point.
        return line.upperPoint
    else:
        x=line.upperPoint.x+(line.lowerPoint.x-line.upperPoint.x)*(point.y-line.upperPoint.y)/(line.lowerPoint.y-line.upperPoint.y)
        return point_2D(x,point.y)

def calSlope(line):
    """Calculate the slope of a given line."""
    deltaX=line.upperPoint.x-line.lowerPoint.x 
    deltaY=line.upperPoint.y-line.lowerPoint.y 

    if deltaX!=0.0:
        return deltaY/deltaX
    else:
        return math.inf

def angleOrder(line1,line2,direction='down'):
    isLess,isGreater,isEqual=False,False,False
    slope1=calSlope(line1)
    slope2=calSlope(line2)

    x1=line1.upperPoint.x-line1.lowerPoint.x
    y1=line1.upperPoint.y-line1.lowerPoint.y
    x2=line2.upperPoint.x-line2.lowerPoint.x
    y2=line2.upperPoint.y-line2.lowerPoint.y

    if abs(crossProductBasic(x1,y1,x2,y2))<=line1.length()*line2.length()*tol:
        # Two lines overlap.
        if line1.halfEdge[0].incidentFace.label<line2.halfEdge[0].incidentFace.label: isLess=True
        else: isGreater=True
        return isLess,isGreater,isEqual
    elif slope1==0.0: isLess=True
    elif slope2==0.0: isGreater=True
    elif slope1*slope2<0 and slope1<slope2: isLess=True # slope1 is negative -- less.
    elif slope1*slope2<0 and slope1>slope2: isGreater=True
    elif slope1<slope2: isGreater=True
    else: isLess=True
    if direction=='down': return isLess,isGreater,isEqual
    elif direction=='up': return isGreater,isLess,isEqual
    else: raise ValueError

def myOrder(line1,line2,point,direction='down'):
    """Return Boolean values: isLess, isGreater, isEqual."""
    if isOnLine(line1,point) and isOnLine(line2,point):
        return angleOrder(line1,line2,direction)

    inter1=calHorizonInter(line1,point)
    inter2=calHorizonInter(line2,point)
    isLess,isGreater,isEqual=False,False,False

    if inter1==inter2:
        x1=line1.upperPoint.x-line1.lowerPoint.x
        y1=line1.upperPoint.y-line1.lowerPoint.y
        x2=line2.upperPoint.x-line2.lowerPoint.x
        y2=line2.upperPoint.y-line2.lowerPoint.y

        if abs(crossProductBasic(x1,y1,x2,y2))<=line1.length()*line2.length()*tol:
            # Two lines overlap.
            if line1.halfEdge[0].incidentFace.label<line2.halfEdge[0].incidentFace.label: isLess=True
            else: isGreater=True
            return isLess,isGreater,isEqual
        elif line1.upperPoint==line2.upperPoint:
            return angleOrder(line1,line2,'down')
        elif line1.lowerPoint==line2.lowerPoint:
            return angleOrder(line1,line2,'up')
        else:
            inter=calIntersection(line1,line2)
            if not inter: 
                print(line1);print(line2)
                inter=calIntersection(line1,line2)
            assert inter,'Error!'
            assert inter<point,'Error!'
            return angleOrder(line1,line2,'up')

    elif inter1<inter2:
        isLess=True
    else:
        isGreater=True

    return isLess,isGreater,isEqual

class element(object):
    # Overlapping element carries not only the element numbering but also the corresponding polygons.
    def __init__(self,theList):
        self.numbering=theList
        self.face=None # The corresponding face in the doubly-connected edge representation.
        self.position=None # Store the position of the element in the input mesh data.
        self.polygons=[]
        self.vertex=[]
    
    def __str__(self):
        return str(self.numbering)

## First step: Find possible elements in different meshes that may overlap.
def elementFilter(coordinates,meshes,transform=False):
    """Here we return the non-overlapping meshes and roughly the possible overlapping meshes.
    
    Input: coordinates and meshes obtained from the input text.
    
    Output: RegularMesh corresponding to non-overlapping elements and the possible overlapping elements in overlappingMesh."""

    numberOfMesh=len(meshes) # Number of meshes.
    overlappingMesh=[None]*numberOfMesh
    squareBound=[None]*numberOfMesh # Ranges of these meshes.

    if transform==False: # The subroutine works as a filter.
        if numberOfMesh==1: return meshes,overlappingMesh # Only one mesh
        
        for i in range(numberOfMesh-1):
            squareBound[i]=findRange(coordinates,meshes[i]) # Find bounds for all meshes except the last one.

        scale[0]=squareBound[0][1]-squareBound[0][0]
        scale[1]=squareBound[0][3]-squareBound[0][2]

        for i in range(numberOfMesh-1,-1,-1):
            xMin,xMax=math.inf,-math.inf
            yMin,yMax=math.inf,-math.inf
            overlapping=[]

            for j in range(len(meshes[i])):
                [xMin1,xMax1,yMin1,yMax1]=findRange(coordinates,[meshes[i][j]]) # Range of this single element

                if any([ifIntersect([xMin1,xMax1,yMin1,yMax1],squareBound[k]) for k in range(numberOfMesh) if k!=i]):
                    # If this element may intersect with any of other meshes.
                    overlapping.append(element(meshes[i][j]))
                    overlapping[-1].position=(i,j)

                    xMin=min(xMin,xMin1)
                    xMax=max(xMax,xMax1)
                    yMin=min(yMin,yMin1)
                    yMax=max(yMax,yMax1) # In this step, we update the range for mesh i, considering only overlapping elements.

                    meshes[i][j]=None # Remove it from the mesh list, since it is already added into the overlapping mesh list.
                else: pass
            
            squareBound[i]=[xMin,xMax,yMin,yMax]
            overlappingMesh[i]=overlapping

    elif transform==True: # All elements are regarded as overlapping elements.
        squareBound=findRange(coordinates,meshes[0])
        scale[0]=squareBound[1]-squareBound[0]
        scale[1]=squareBound[3]-squareBound[2]

        for i in range(numberOfMesh):
            overlapping=[]

            for j in range(len(meshes[i])):
                overlapping.append(element(meshes[i][j]))
                overlapping[-1].position=(i,j)

                meshes[i][j]=None
            
            overlappingMesh[i]=overlapping
            
    else: raise ValueError
    
    return meshes,overlappingMesh

def ifIntersect(range1,range2):
    """Given two squares, test if they intersect.

    Input: Two squares presented by range1 and range2, both formated as [x_min, x_max, y_min, y_max]."""

    if range1[1]<=range2[0] or range1[0]>=range2[1] or range1[2]>=range2[3] or range1[3]<=range2[2]: 
        return False
    else: return True

def findRange(coordinates,mesh):
    """Calculate the square range of a mesh.
    
    Output: [x_min, x_max, y_min, y_max]."""

    if len(mesh)>1:
        (xMin1,xMax1,yMin1,yMax1)=findRange(coordinates,mesh[0 : len(mesh)//2])
        (xMin2,xMax2,yMin2,yMax2)=findRange(coordinates,mesh[len(mesh)//2 : len(mesh)])
        xMin=min(xMin1,xMin2)
        xMax=max(xMax1,xMax2)
        yMin=min(yMin1,yMin2)
        yMax=max(yMax1,yMax2)

    else:
        xMax=xMin=coordinates[mesh[0][0]][0]
        yMax=yMin=coordinates[mesh[0][0]][1]

        if len(mesh[0])==4 or len(mesh[0])==9: numberOfNode=4 # Quadrilateral
        elif len(mesh[0])==3 or len(mesh[0])==6: numberOfNode=3 # Trianle
        else: raise ValueError

        for i in range(1,numberOfNode):
            xMax=max(xMax,coordinates[mesh[0][i]][0])
            xMin=min(xMin,coordinates[mesh[0][i]][0])
            yMax=max(yMax,coordinates[mesh[0][i]][1])
            yMin=min(yMin,coordinates[mesh[0][i]][1])

    return [xMin,xMax,yMin,yMax]

def insertEventQueue(node,k):
    """Insert a point into the event queue. If the coordinates have already been added, merge the point with the existing point."""

    # This function is different from the one in bstAVL because we do not want to have repeated event points.
    if node.key: # node.key is None only when the tree is empty.
        if k==node.key:
            # The coordinates have already been added into the BST.
            # Extend the incidentL and incidentU lists.
            # Merge vertices that have same coordinates, to keep track of involvedMesh and leftList.
            if k.num is not None: # k.num is not None only when the point is generated from toDC.toDoublyConnected.
                # Then, k.label is not None.
                node.key.overlappingVertex[k.label]=k

            # Since k is not used in establishing the mesh overlay, the coordinates are just floats.
            k.x=float(k.x)
            k.y=float(k.y)

            node.key.incidentL.extend(k.incidentL)
            for i in k.incidentL: 
                i.lowerPoint=node.key
                if i.halfEdge[0].origin==node.key: 
                    i.halfEdge[0].origin=node.key
                    i.halfEdge[1].destination=node.key
                else:
                    i.halfEdge[0].destination=node.key
                    i.halfEdge[1].origin=node.key

            node.key.incidentU.extend(k.incidentU)
            for i in k.incidentU:
                i.upperPoint=node.key
                if i.halfEdge[0].origin==node.key: 
                    i.halfEdge[0].origin=node.key
                    i.halfEdge[1].destination=node.key
                else:
                    i.halfEdge[0].destination=node.key
                    i.halfEdge[1].origin=node.key

            if len(node.key.incidentC)==0 and len(k.incidentC)!=0: node.key.incidentC.append(k.incidentC[0]) 
                # If node.key is an intersection, incidentC will not be empty.
        elif k>node.key:
            if node.right:
                insertEventQueue(node.right,k)
            else:
                node.right=bstAVL.new_node(k)
                node.right.parent=node
                bstAVL.maintain(node)
        elif k<node.key:
            if node.left:
                insertEventQueue(node.left,k)
            else:
                node.left=bstAVL.new_node(k)
                node.left.parent=node
                bstAVL.maintain(node)
    else:
        node.key=k

def meshOverlay(lineList,vertexList,halfEdgeList,faceList,numberOfMesh):
    """Calculate the new subdivision given the overlapping meshes and their doubly-connected edge representation."""
    # Initiate the polygon list.
    polygons=[]

    # Initiate the event queue.
    eventQueue=bstAVL.new_node(None)

    for i in range(len(vertexList)):
        insertEventQueue(eventQueue,vertexList[i])
    
    # Initiate the status (The total status).
    status=bstAVLLines.new_node(None)

    # We need some more status trees, one for each mesh.
    statusList=[None]*numberOfMesh
    for i in range(numberOfMesh):
        statusList[i]=bstAVLLines.new_node(None)

    # We store the list of event points.
    eventList=[]

    while eventQueue.key: # Do while the event queue is not empty.
        nextEvent=bstAVL.find_max(eventQueue)
        point=nextEvent.key
        eventList.append(point)
        bstAVL.delete(nextEvent) # Remove it from the event queue.

        point.involvedMesh=[0]*numberOfMesh # We need the information how many meshes are involved.
        point.leftList=[None]*numberOfMesh # This will be used to store information about what elements the new faces belong to.

        # At each event point, we need to do some work.
        # print("The event point is %s."%point)
        handleEventPointInMesh(point,status,statusList,halfEdgeList,lineList,eventQueue,numberOfMesh)

    # Now we merge meshes at each event point.
    # This step is not finished in handleEventPointInMesh to avoid round-off errors, e.g. crossProduct==0...
    # Because, after handleEventPointInMesh, now overlapping edges must have exactly the same origin and destination.
    for i in eventList:
        buildConnectivity(i)
        # pass

    # Assemble halfEdges to obtain new faces.
    # Travese to calculate all new faces.
    boundary=[]

    for i in halfEdgeList:
        if i.newFace is None:
            # Create a new face representing this polygon.
            newFace=toDC.face(None,i)
            i.newFace=newFace
            newFace.incidentElements=[None]*numberOfMesh
            # Find all elements that contain this polygon.
            for j in i.faceList: 
                if j:
                    # Meshes that contain this sub-region:
                    newFace.incidentElements[j.label]=j.element
                    if not isZeroArea(newFace):
                        j.element.polygons.append(newFace)
            
            for j in range(numberOfMesh):
                if i.origin.involvedMesh[j]==0 and i.origin.leftList[j]:
                    if i.origin.leftList[j].incidentFace: 
                        # Meshes that contain this sub-region:
                        newFace.incidentElements[j]=i.origin.leftList[j].incidentFace.element
                        if not isZeroArea(newFace):
                            i.origin.leftList[j].incidentFace.element.polygons.append(newFace)
            
            # Find all half-edges incident to this new face.
            nextEdge=i.next
            while nextEdge.newFace is None:
                nextEdge.newFace=newFace
                nextEdge=nextEdge.next
            
            # If the new face is not the exterior space and the new face is not formed by two half-edges.
            if (not isOuterLoop(newFace)) and (not isZeroArea(newFace)):
                polygons.append(newFace)
            elif isOuterLoop(newFace): boundary.append(newFace.outer)

    # Find mesh boundaries for constructing weight functions.
    for i in lineList:
        # Convert the coordinates to floats.
        toFloatLine(i)
        # If the two sides of this line belong to different elements, it is a boundary of some mesh.
        leftFace=findFace(i.halfEdge[0])
        rightFace=findFace(i.halfEdge[1])

        if (leftFace) and (rightFace):
            for j in range(numberOfMesh):
                if ((leftFace.incidentElements[j] is None) and rightFace.incidentElements[j]) or \
                    (leftFace.incidentElements[j] and (rightFace.incidentElements[j] is None)):
                    if i.upperPoint.overlappingVertex[j]: i.upperPoint.overlappingVertex[j].weight=0
                    if i.lowerPoint.overlappingVertex[j]: i.lowerPoint.overlappingVertex[j].weight=0

    return polygons, boundary

def toFloatLine(line):
    """Convert the coordinates of vertices of a line segment to floats."""
    line.upperPoint.x=float(line.upperPoint.x)
    line.upperPoint.y=float(line.upperPoint.y)
    line.lowerPoint.x=float(line.lowerPoint.x)
    line.lowerPoint.y=float(line.lowerPoint.y)

def isOuterLoop(newFace):
    """NewFace is in the exterior iff all its incident elements are None."""

    for i in newFace.incidentElements:
        if i: return False
    
    return True

def isZeroArea(newFace):
    """If newFace has only two incident edges, it can be removed."""

    if newFace.outer.next.next is newFace.outer: return True
    else: return False

def findFace(halfEdge):
    """Find the first face on the left of halfEdge that is not degenerated."""

    while halfEdge.next.next is halfEdge:
        halfEdge=halfEdge.next.twin

    for i in halfEdge.newFace.incidentElements:
        if i: return halfEdge.newFace
    
    return None

def edgeSplit(point,edge,halfEdgeList,lineList):
    """Split a line segment at its interior event point."""
    # The event point must be an interior point of the edge.
    oldUpperPoint=edge.upperPoint

    edge.upperPoint=point
    point.incidentU.append(edge) # We don't modify point.incidentL because it will no longer be used.

    # The half-edges incident to the (new) upper end point should remain unchanged.
    if edge.halfEdge[0].origin is edge.lowerPoint:
        edge.halfEdge[0].origin=point
        edge.halfEdge[1].destination=point

        # Create two new half edges and add them to the half-edge list.
        newHalfEdge=toDC.halfEdge(edge.lowerPoint,edge.halfEdge[0].incidentFace) # halfEdge[0] must have an incident face.
        newHalfEdge.destination=point
        otherNewHalfEdge=toDC.halfEdge(point,edge.halfEdge[1].incidentFace)
        otherNewHalfEdge.destination=edge.lowerPoint

        # Create a new line.
        newLine=line(oldUpperPoint,point)
        newLine.halfEdge=edge.halfEdge
        edge.halfEdge=[newHalfEdge,otherNewHalfEdge]
        newHalfEdge.line=edge
        otherNewHalfEdge.line=edge
        newLine.halfEdge[0].line=newLine
        newLine.halfEdge[1].line=newLine

        if point.incidentEdge is None: point.incidentEdge=newLine.halfEdge[0]

        newHalfEdge.prev=newLine.halfEdge[0].prev
        newLine.halfEdge[0].prev.next=newHalfEdge
        newHalfEdge.next=newLine.halfEdge[0]
        newLine.halfEdge[0].prev=newHalfEdge

        otherNewHalfEdge.prev=newLine.halfEdge[1]
        otherNewHalfEdge.next=newLine.halfEdge[1].next
        newLine.halfEdge[1].next.prev=otherNewHalfEdge
        newLine.halfEdge[1].next=otherNewHalfEdge

    else:
        edge.halfEdge[0].destination=point
        edge.halfEdge[1].origin=point

        # Create two new half edges and add them to the half-edge list.
        newHalfEdge=toDC.halfEdge(point,edge.halfEdge[0].incidentFace)
        newHalfEdge.destination=edge.lowerPoint
        otherNewHalfEdge=toDC.halfEdge(edge.lowerPoint,edge.halfEdge[1].incidentFace)
        otherNewHalfEdge.destination=point

        newLine=line(oldUpperPoint,point)
        newLine.halfEdge=edge.halfEdge
        edge.halfEdge=[newHalfEdge,otherNewHalfEdge]
        newHalfEdge.line=edge
        otherNewHalfEdge.line=edge
        newLine.halfEdge[0].line=newLine
        newLine.halfEdge[1].line=newLine

        if point.incidentEdge is None: point.incidentEdge=newLine.halfEdge[1]

        newHalfEdge.prev=newLine.halfEdge[0]
        newHalfEdge.next=newLine.halfEdge[0].next
        newLine.halfEdge[0].next.prev=newHalfEdge
        newLine.halfEdge[0].next=newHalfEdge

        otherNewHalfEdge.next=newLine.halfEdge[1]
        otherNewHalfEdge.prev=newLine.halfEdge[1].prev
        newLine.halfEdge[1].prev.next=otherNewHalfEdge
        newLine.halfEdge[1].prev=otherNewHalfEdge

    halfEdgeList.append(newHalfEdge)
    halfEdgeList.append(otherNewHalfEdge)    
    lineList.append(newLine)
    newHalfEdge.twin=otherNewHalfEdge
    otherNewHalfEdge.twin=newHalfEdge

def handleEventPointInMesh(point,status,statusList,halfEdgeList,lineList,eventQueue,numberOfMesh):
    """Handle each event point."""
    # First step: Find all lines that contain the event point.
    # We can start from the first element and find the remaining all neighborhoods.
    if len(point.incidentL)+len(point.incidentC)==0: 
        findNeighbor=None
    else: 
        if point.incidentL: 
            startingLine=point.incidentL[0]
            point.incidentC=[]
            # Since this is an end point, we do not split the edge.
        else: 
            startingLine=point.incidentC[0]
            point.incidentC=[] # We no longer need incidentC, all incidentC edges will be added into incidentU.

            # We split the edge into two.
            edgeSplit(point,startingLine,halfEdgeList,lineList)
        
        findNeighbor=bstAVLLines.find(status,startingLine,point,'up')
        # Find the position of the starting line in the BST.
        findAgain=bstAVLLines.find(statusList[startingLine.halfEdge[0].incidentFace.label],startingLine,point,'up')
        # Also find its position in the corresponding sub-BST.
        assert (findNeighbor and findAgain),"Some line not found!"
        # The first incident edge is not on the outer boundary, so it will always have a label.
        point.involvedMesh[startingLine.halfEdge[0].incidentFace.label]=1
        
        if startingLine.halfEdge[0].origin is point:
            point.leftList[startingLine.halfEdge[0].incidentFace.label]=startingLine.halfEdge[0]
        else: point.leftList[startingLine.halfEdge[0].incidentFace.label]=startingLine.halfEdge[1]
        # The edges in leftList will always have point as the common origin if they are incident edges.

    leftKey=rightKey=None

    if findNeighbor: # Delete all incidentL and incidentC from the BST, and find all incidentC lines (split them and add into incidentU).
        left=bstAVLLines.predecessor(findNeighbor)
        right=bstAVLLines.successor(findNeighbor)

        if left: leftKey=left.key 
        if right: rightKey=right.key 

        bstAVLLines.delete(findNeighbor)
        bstAVLLines.delete(findAgain) # Also delete it from the sub-BST.

        left=None
        if leftKey: 
            left=bstAVLLines.find(status,leftKey,point,'up')
            anotherLeft=bstAVLLines.find(statusList[leftKey.halfEdge[0].incidentFace.label],leftKey,point,'up')
            assert (left and anotherLeft),"Some line not found!"

        while left and isOnLine(left.key,point):
            if left.key.upperPoint!=point and left.key.lowerPoint!=point: # This is an incidentC edge.
                edgeSplit(point,leftKey,halfEdgeList,lineList)

            point.involvedMesh[leftKey.halfEdge[0].incidentFace.label]=1
            if point.leftList[leftKey.halfEdge[0].incidentFace.label] is None:
                if leftKey.halfEdge[0].origin is point:
                    point.leftList[leftKey.halfEdge[0].incidentFace.label]=leftKey.halfEdge[0]
                else: point.leftList[leftKey.halfEdge[0].incidentFace.label]=leftKey.halfEdge[1]
            
            newLeft=bstAVLLines.predecessor(left)
            
            newLeftKey=None
            if newLeft: newLeftKey=newLeft.key
            
            bstAVLLines.delete(left)
            bstAVLLines.delete(anotherLeft)

            left=None
            if newLeftKey:
                left=bstAVLLines.find(status,newLeftKey,point,'up')
                anotherLeft=bstAVLLines.find(statusList[newLeftKey.halfEdge[0].incidentFace.label],newLeftKey,point,'up')
                assert (left and anotherLeft),"Some line not found!"
            leftKey=newLeftKey

        right=None
        if rightKey: 
            right=bstAVLLines.find(status,rightKey,point,'up')
            anotherRight=bstAVLLines.find(statusList[rightKey.halfEdge[0].incidentFace.label],rightKey,point,'up')
            assert (right and anotherRight),"Some line not found!"

        while right and isOnLine(right.key,point):
            if right.key.upperPoint!=point and right.key.lowerPoint!=point:
                edgeSplit(point,rightKey,halfEdgeList,lineList)

            point.involvedMesh[rightKey.halfEdge[0].incidentFace.label]=1
            if point.leftList[rightKey.halfEdge[0].incidentFace.label] is None:
                if rightKey.halfEdge[0].origin is point:
                    point.leftList[rightKey.halfEdge[0].incidentFace.label]=rightKey.halfEdge[0]
                else: point.leftList[rightKey.halfEdge[0].incidentFace.label]=rightKey.halfEdge[1]

            newRight=bstAVLLines.successor(right)

            newRightKey=None
            if newRight: newRightKey=newRight.key 

            bstAVLLines.delete(right)
            bstAVLLines.delete(anotherRight)

            right=None
            if newRightKey: 
                right=bstAVLLines.find(status,newRightKey,point,'up')
                anotherRight=bstAVLLines.find(statusList[newRightKey.halfEdge[0].incidentFace.label],newRightKey,point,'up')
                assert (right and anotherRight),"Some line not found!"
            rightKey=newRightKey
    
    for i in point.incidentU:
        bstAVLLines.insert(status,i,point)

        # bstAVLLines.checkBST(status,point)
        # Also update the involvedMesh list and leftList.
        point.involvedMesh[i.halfEdge[0].incidentFace.label]=1
        if point.leftList[i.halfEdge[0].incidentFace.label] is None:
            if i.halfEdge[0].origin is point:
                point.leftList[i.halfEdge[0].incidentFace.label]=i.halfEdge[0]
            else: point.leftList[i.halfEdge[0].incidentFace.label]=i.halfEdge[1]

        # These edges are also added into those sub-BST's.
        bstAVLLines.insert(statusList[i.halfEdge[0].incidentFace.label],i,point)

    if len(point.incidentU)==0: # This event point is not an upper end point for any line segment.
        left=None
        right=None

        if leftKey: left=bstAVLLines.find(status,leftKey,point)
        if rightKey: right=bstAVLLines.find(status,rightKey,point)
        
        findNewEvent(left,right,point,eventQueue) # Calculate the intersection if it exists.
    else:
        startingLine=point.incidentU[0]

        leftMost=rightMost=bstAVLLines.find(status,startingLine,point)

        while True:
            newRight=None
            if rightMost: newRight=bstAVLLines.successor(rightMost)
            if newRight and isOnLine(newRight.key,point):
                if point!=newRight.key.upperPoint and point!=newRight.key.lowerPoint:
                    # This is the hidden incidentC line segment!
                    edgeSplit(point,newRight.key,halfEdgeList,lineList)
                    point.involvedMesh[newRight.key.halfEdge[0].incidentFace.label]=1
                    if point.leftList[newRight.key.halfEdge[0].incidentFace.label] is None:
                        if newRight.key.halfEdge[0].origin is point:
                            point.leftList[newRight.key.halfEdge[0].incidentFace.label]=newRight.key.halfEdge[0]
                        else: point.leftList[newRight.key.halfEdge[0].incidentFace.label]=newRight.key.halfEdge[1]
                rightMost=newRight
            else:
                break
        while True:
            newLeft=None
            if leftMost: newLeft=bstAVLLines.predecessor(leftMost)
            if newLeft and isOnLine(newLeft.key,point):
                if point!=newLeft.key.upperPoint and point!=newLeft.key.lowerPoint:
                    edgeSplit(point,newLeft.key,halfEdgeList,lineList)
                    point.involvedMesh[newLeft.key.halfEdge[0].incidentFace.label]=1
                    if point.leftList[newLeft.key.halfEdge[0].incidentFace.label] is None:
                        if newLeft.key.halfEdge[0].origin is point:
                            point.leftList[newLeft.key.halfEdge[0].incidentFace.label]=newLeft.key.halfEdge[0]
                        else: point.leftList[newLeft.key.halfEdge[0].incidentFace.label]=newLeft.key.halfEdge[1]
                leftMost=newLeft
            else:
                break
        
        findNewEvent(newLeft,leftMost,point,eventQueue)
        findNewEvent(rightMost,newRight,point,eventQueue)

    # For meshes that are not involved at this event point, we need also to store the immediately left edge. 
    tempNode=point_2D(point.x,point.y+1.0)
    # Create a new edge that is parallel to the y-axis.
    tempEdge=line(point,tempNode)

    for i in range(numberOfMesh):
        if point.leftList[i] is None: # This mesh is not involved.
            bstAVLLines.insert(statusList[i],tempEdge,point)
            findTemp=bstAVLLines.find(statusList[i],tempEdge,point)
            immediateLeft=bstAVLLines.successor(findTemp)

            if immediateLeft: 
                if immediateLeft.key.halfEdge[0].origin>immediateLeft.key.halfEdge[0].destination:
                    point.leftList[i]=immediateLeft.key.halfEdge[0]
                else: point.leftList[i]=immediateLeft.key.halfEdge[1]
            bstAVLLines.delete(findTemp) # Delete this temp edge from the status.

def buildConnectivity(point):
    """Correctly connect half-edges around each event point."""
    # Now we need to correctly establish the connections between half edges at the event point.

    for i in range(len(point.involvedMesh)): 
        if point.involvedMesh[i]: 
            startingNumber=i
            break
    
    # An edge in the first involved mesh.
    someEdge=point.leftList[startingNumber]

    while True:
        if len(someEdge.faceList)==0:
            someEdge.faceList.append(someEdge.incidentFace)
            someEdge=someEdge.twin.next
        else: break
    
    if sum(point.involvedMesh)>1: # At least two meshes are involved.
        for i in range(startingNumber+1,len(point.involvedMesh)):
            if point.involvedMesh[i]: mergeMesh(point.leftList[startingNumber],point.leftList[i])

def isLeftOverlapping(halfEdge1,halfEdge2):
    """Given two overlapping half edges (having same origin), we further define their relative position."""

    # halfEdge1 and halfEdge2 must have different labels, since edges in the same mesh cannot overlap.
    if halfEdge1.origin>halfEdge1.destination: # Origin is an upper end point.
        if getLabel(halfEdge1)<getLabel(halfEdge2): return False # We say halfEdge1 is on the left of halfEdge2.
        else: return True
    else: 
        if getLabel(halfEdge1)<getLabel(halfEdge2): return True
        else: return False

def crossProduct(v1,v2):
    """Given two half-edges that have common origin, calculate their cross product."""
    x1=v1.destination.x-v1.origin.x
    y1=v1.destination.y-v1.origin.y 

    x2=v2.destination.x-v2.origin.x 
    y2=v2.destination.y-v2.origin.y 

    return crossProductBasic(x1,y1,x2,y2)

def crossProductBasic(x1,y1,x2,y2):
    return x1*y2-x2*y1

def dotProduct(v1,v2):
    """Given two half-edges that have common origin, calculate their dot product."""
    x1=v1.destination.x-v1.origin.x
    y1=v1.destination.y-v1.origin.y 

    x2=v2.destination.x-v2.origin.x 
    y2=v2.destination.y-v2.origin.y 

    return dotProductBasic(x1,y1,x2,y2)

def dotProductBasic(x1,y1,x2,y2):
    return x1*x2+y1*y2

def getLabel(halfEdge):
    """Return the mesh number of halfEdge."""

    if halfEdge.incidentFace:
        return halfEdge.incidentFace.label 
    else:
        return halfEdge.twin.incidentFace.label 

def isWellOrdered(halfEdge1,halfEdge2,halfEdge3):
    """Given three half-edges that have the same origin, check if they are in good cyclic order: halfEdge1<halfEdge2<halfEdge3."""
    # Here we use "==0" because either two edges totally overlap or not co-linear.

    if crossProduct(halfEdge1,halfEdge3)<0: # The angle between these two edges is less than pi.
        if crossProduct(halfEdge1,halfEdge2)<0 and crossProduct(halfEdge2,halfEdge3)<0:
            return True
        elif crossProduct(halfEdge1,halfEdge2)==0 and dotProduct(halfEdge1,halfEdge2)>0:
            if isLeftOverlapping(halfEdge1,halfEdge2): return True
        elif crossProduct(halfEdge2,halfEdge3)==0 and dotProduct(halfEdge2,halfEdge3)>0:
            if isLeftOverlapping(halfEdge2,halfEdge3): return True

    elif crossProduct(halfEdge1,halfEdge3)>0: # The angle between these two edges is greater than pi.
        if crossProduct(halfEdge1,halfEdge2)<0 or crossProduct(halfEdge2,halfEdge3)<0:
            return True
        elif crossProduct(halfEdge1,halfEdge2)==0 and dotProduct(halfEdge1,halfEdge2)>0:
            if isLeftOverlapping(halfEdge1,halfEdge2): return True
        elif crossProduct(halfEdge1,halfEdge2)==0 and dotProduct(halfEdge1,halfEdge2)<0:
            return True
        elif crossProduct(halfEdge2,halfEdge3)==0 and dotProduct(halfEdge2,halfEdge3)>0:
            if isLeftOverlapping(halfEdge2,halfEdge3): return True
        elif crossProduct(halfEdge2,halfEdge3)==0 and dotProduct(halfEdge2,halfEdge3)<0:
            return True

    elif crossProduct(halfEdge1,halfEdge3)==0 and dotProduct(halfEdge1,halfEdge3)<0: # The angle is pi.
        if crossProduct(halfEdge1,halfEdge2)<0: return True
        elif crossProduct(halfEdge1,halfEdge2)==0:
            if dotProduct(halfEdge1,halfEdge2)>0: 
                if isLeftOverlapping(halfEdge1,halfEdge2): return True
            else:
                if isLeftOverlapping(halfEdge2,halfEdge3): return True
    
    return False

def mergeMesh(halfEdge1,halfEdge2):
    """halfEdge1 gives the chain of already sorted half-edges;
    halfEdge2 gives a new chain to be merged into the existing list.
    We require that the label number of halfEdge2 is larger."""

    count=0
    # First, find the cyclic position of halfEdge1 in the chain of halfEdge2.
    while not isWellOrdered(halfEdge2,halfEdge1,halfEdge2.twin.next):
        halfEdge2=halfEdge2.twin.next
        count+=1
        if count>100: raise ValueError
    
    # Find the starting edge.
    while isWellOrdered(halfEdge2,halfEdge1.prev.twin,halfEdge2.twin.next) and \
        isWellOrdered(halfEdge1.prev.twin,halfEdge1,halfEdge2.twin.next):
        halfEdge1=halfEdge1.prev.twin
        count+=1
        if count>100: raise ValueError
    
    startingEdge=halfEdge1
    startingEdge2=halfEdge2.twin.next

    nextHalfEdge2=halfEdge2.twin.next

    halfEdge1=halfEdge1.twin.next
    nextHalfEdge1=halfEdge1.twin.next

    while halfEdge1 is not startingEdge:
        while isWellOrdered(halfEdge2,halfEdge1,nextHalfEdge2):
            # If the next half-edge enter the same interval from the other direction.
            if isWellOrdered(halfEdge2,nextHalfEdge1,halfEdge1):
                halfEdge1.faceList.append(nextHalfEdge2.incidentFace)
                halfEdge1=nextHalfEdge1
                nextHalfEdge1=halfEdge1.twin.next

                tempEdge=nextHalfEdge2

                while tempEdge is not halfEdge2:
                    tempEdge.faceList.extend(halfEdge1.faceList)
                    tempEdge.faceList.append(tempEdge.incidentFace)
                    tempEdge=tempEdge.twin.next

                    count+=1
                    if count>100: raise ValueError

                halfEdge2.faceList.extend(halfEdge1.faceList)
                halfEdge2.faceList.append(halfEdge2.incidentFace)

                break

            halfEdge1.faceList.append(nextHalfEdge2.incidentFace)
            halfEdge1=nextHalfEdge1
            nextHalfEdge1=halfEdge1.twin.next

            count+=1
            if count>100: raise ValueError
        
        halfEdge1.prev.next=nextHalfEdge2
        nextHalfEdge2.prev=halfEdge1.prev

        while not isWellOrdered(halfEdge2,halfEdge1,nextHalfEdge2):
            halfEdge2=nextHalfEdge2
            nextHalfEdge2=halfEdge2.twin.next
            halfEdge2.faceList.extend(halfEdge1.faceList)
            halfEdge2.faceList.append(halfEdge2.incidentFace)

            count+=1
            if count>100: raise ValueError

        halfEdge2.twin.next=halfEdge1
        halfEdge1.prev=halfEdge2.twin

    startingEdge.faceList.append(startingEdge2.incidentFace)

def traversalEdges(halfEdge):
    """For debugging."""
    print(halfEdge)
    startingEdge=halfEdge
    halfEdge=halfEdge.twin.next

    while halfEdge is not startingEdge:
        print(halfEdge)
        halfEdge=halfEdge.twin.next

def calIntersection(lineLeft,lineRight):
    """Calculate the intersection of two line segments if it exists."""
    x1=lineLeft.upperPoint.x-lineLeft.lowerPoint.x
    y1=lineLeft.upperPoint.y-lineLeft.lowerPoint.y
    x2=lineRight.upperPoint.x-lineRight.lowerPoint.x
    y2=lineRight.upperPoint.y-lineRight.lowerPoint.y

    if abs(crossProductBasic(x1,y1,x2,y2))<=lineLeft.length()*lineRight.length()*tol:
        return None
        # If two edges overlap, we do not need to calculate the intersection.
    elif lineLeft.lowerPoint==lineRight.lowerPoint:
        return None
    else:
        Px=((lineLeft.upperPoint.x*lineLeft.lowerPoint.y-lineLeft.upperPoint.y*lineLeft.lowerPoint.x)\
            *(lineRight.upperPoint.x-lineRight.lowerPoint.x)-(lineLeft.upperPoint.x-lineLeft.lowerPoint.x)\
                *(lineRight.upperPoint.x*lineRight.lowerPoint.y-lineRight.upperPoint.y*lineRight.lowerPoint.x))\
                    /((lineLeft.upperPoint.x-lineLeft.lowerPoint.x)*(lineRight.upperPoint.y-lineRight.lowerPoint.y)-\
                        (lineLeft.upperPoint.y-lineLeft.lowerPoint.y)*(lineRight.upperPoint.x-lineRight.lowerPoint.x))
        Py=((lineLeft.upperPoint.x*lineLeft.lowerPoint.y-lineLeft.upperPoint.y*lineLeft.lowerPoint.x)\
            *(lineRight.upperPoint.y-lineRight.lowerPoint.y)-(lineLeft.upperPoint.y-lineLeft.lowerPoint.y)\
                *(lineRight.upperPoint.x*lineRight.lowerPoint.y-lineRight.upperPoint.y*lineRight.lowerPoint.x))\
                    /((lineLeft.upperPoint.x-lineLeft.lowerPoint.x)*(lineRight.upperPoint.y-lineRight.lowerPoint.y)-\
                        (lineLeft.upperPoint.y-lineLeft.lowerPoint.y)*(lineRight.upperPoint.x-lineRight.lowerPoint.x))
        
        # Dealing with round-off errors.
        if (lineLeft.upperPoint.y-lineLeft.lowerPoint.y)<scale[1]*tol: Py=lineLeft.upperPoint.y
        if (lineRight.upperPoint.y-lineRight.lowerPoint.y)<scale[1]*tol: Py=lineRight.upperPoint.y
        if abs(lineLeft.upperPoint.x-lineLeft.lowerPoint.x)<scale[0]*tol: Px=lineLeft.upperPoint.x
        if abs(lineRight.upperPoint.x-lineRight.lowerPoint.x)<scale[0]*tol: Px=lineRight.upperPoint.x

        intersection=point_2D(Px,Py)
        if intersection==lineLeft.upperPoint: 
            intersection.x=lineLeft.upperPoint.x
            intersection.y=lineLeft.upperPoint.y
        if intersection==lineLeft.lowerPoint: 
            intersection.x=lineLeft.lowerPoint.x
            intersection.y=lineLeft.lowerPoint.y
        if intersection==lineRight.upperPoint: 
            intersection.x=lineRight.upperPoint.x
            intersection.y=lineRight.upperPoint.y
        if intersection==lineRight.lowerPoint: 
            intersection.x=lineRight.lowerPoint.x
            intersection.y=lineRight.lowerPoint.y
        
        if lineLeft.upperPoint.y==lineLeft.lowerPoint.y: intersection.y=lineLeft.upperPoint.y
        if lineRight.upperPoint.y==lineRight.lowerPoint.y: intersection.y=lineRight.upperPoint.y

        if isOnLine(lineLeft,intersection) and isOnLine(lineRight,intersection):
            if lineLeft.upperPoint!=intersection and lineLeft.lowerPoint!=intersection:
                intersection.incidentC.append(lineLeft)
            if lineRight.upperPoint!=intersection and lineRight.lowerPoint!=intersection:
                intersection.incidentC.append(lineRight)
            return intersection
        else: return None

def findNewEvent(left,right,point,eventQueue):
    """Given two line segments, return their intersection if it exists. And add it into the event queue."""
    if left and right:
        lineLeft=left.key
        lineRight=right.key
        intersection=calIntersection(lineLeft,lineRight)

        if intersection and intersection<point:
            intersection.overlappingVertex=[None]*len(point.involvedMesh)
            insertEventQueue(eventQueue,intersection)

def rhoCalculation(polygons):
    # Now we also add weights for the meshes whose numbers != 0.
    for i in polygons:
        if triangulation.notOverlap(i): continue

        startingEdge=i.outer
        edge=i.outer

        while True:
            if edge.origin.rho is None:
                edge.origin.rho=np.zeros(len(i.incidentElements))
                xy=np.array([edge.origin.x,edge.origin.y])
                
                for j in i.incidentElements:
                    if j:
                        coord=[]
                        weight=[]
                        for k in j.vertex:
                            coord.append([k.x,k.y])
                            weight.append(k.weight)
                        coord=np.array(coord)
                        weight=np.array(weight)

                        if len(j.numbering)==4 or len(j.numbering)==9:
                            isoCoord=invShapeFunction.invQuad(coord,xy)
                            Nmat,_=invShapeFunction.quadShapeFunction(isoCoord)
                            edge.origin.rho[j.position[0]]=np.inner(Nmat.reshape(-1),weight)
                        elif len(j.numbering)==3 or len(j.numbering)==6:
                            isoCoord=invShapeFunction.invTri(coord,xy)
                            edge.origin.rho[j.position[0]]=np.inner(isoCoord,weight)
                        else: raise ValueError
                
                # Now if we have two meshes completely overlapping, then the second mesh would have rho=0.9.
                for j in range(1,len(edge.origin.rho)):
                    edge.origin.rho[j] *= 9 # It should be 9.

                edge.origin.rho/=sum(edge.origin.rho)

            edge=edge.next
            if edge is startingEdge: break

if __name__=='__main__':
    pass