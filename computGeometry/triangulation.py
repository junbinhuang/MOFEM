import math
import sys, os
sys.path.append(os.path.dirname(sys.path[0]))
from meshTools import toDC

def polygonTriangulation(face):
    """Given a polygonal face, calculate its triangulation.

    Output: If this face is to be triangulated, create a half-edge list, 
    in which each half edge represents a triangle, and attach the list to
    face.triangles."""
    numberOfEdges=countEdges(face.outer)

    # If this element remains intact, we do not triangulate it.
    if notOverlap(face):
        # This element does not overlap with other meshes.
        pass
    else:
        # Triangulate it!
        # Because the number of edges is usually very small (<=6), we use the naive triangulation algorithm.
        # For practical cases, it is necessary to implement also the O(nlogn) triangulation algorithm.
        face.triangles=simpleTriangulation(face.outer,numberOfEdges)

def notOverlap(face):
    # Count the number of incident elements. Also find the first incident element.
    numberOfElements=0
    firstIncidentElement=None
    for i in face.incidentElements:
        if i:
            numberOfElements+=1
            if firstIncidentElement is None: firstIncidentElement=i  
    
    if numberOfElements==1 and len(firstIncidentElement.polygons)==1:
        return True
    else: return False

def countEdges(halfEdge):
    """Given a half-edge of a polygon, count the number of edges in this polygon."""
    numberOfEdges=1
    startingEdge=halfEdge
    edge=startingEdge.next

    while edge is not startingEdge:
        edge=edge.next
        numberOfEdges+=1 

    return numberOfEdges 

def simpleTriangulation(halfEdge,numberOfEdges=None):
    """The naive triangulation for simple polygons.
    It works best when the number of vertices is not too large.
    The cost for the worst case is O(n^2), where n is the number of edges."""
    
    if numberOfEdges is None:
        numberOfEdges=countEdges(halfEdge)

    if numberOfEdges==3:
        return [halfEdge]
    
    # Find the up-left-most vertex.
    upLeftEdge=halfEdge

    for i in range(numberOfEdges-1):
        halfEdge=halfEdge.next
        if halfEdge.origin>upLeftEdge.origin: upLeftEdge=halfEdge
    
    # Check if this diagonal can be drawn.
    interiorPoints=[]
    distance=[]

    triangle=[upLeftEdge.origin,upLeftEdge.next.origin,upLeftEdge.prev.origin,upLeftEdge.origin]

    testEdge=upLeftEdge.next.next

    while testEdge is not upLeftEdge.prev:
        if isInPolygon(testEdge.origin,triangle):
            interiorPoints.append(testEdge)
            distance.append(isLeft(upLeftEdge.next.origin,upLeftEdge.prev.origin,testEdge.origin))
        testEdge=testEdge.next

    # Split the polygon into two polygons. And calculate their triangulations.
    if len(interiorPoints)==0:
        edge1=upLeftEdge.prev
        edge2=upLeftEdge.next

        newHalfEdge=toDC.halfEdge(edge2.origin,None)
        newHalfEdge.destination=edge1.origin
        newHalfEdge.prev=edge2.prev
        newHalfEdge.next=edge1

        otherHalfEdge=toDC.halfEdge(edge1.origin,None)
        otherHalfEdge.destination=edge2.origin
        otherHalfEdge.prev=edge1.prev
        otherHalfEdge.next=edge2

        newHalfEdge.twin=otherHalfEdge
        otherHalfEdge.twin=newHalfEdge

        edge1.prev.next=otherHalfEdge
        edge1.prev=newHalfEdge

        edge2.prev.next=newHalfEdge
        edge2.prev=otherHalfEdge

        triangleList1=simpleTriangulation(edge1,3)
        triangleList2=simpleTriangulation(edge2,numberOfEdges-1)

    else:
        edge1=upLeftEdge

        maxDistance=distance[0]
        maxIndex=0
        for i in range(len(distance)):
            if distance[i]>maxDistance: 
                maxDistance=distance[i]
                maxIndex=i
        
        edge2=interiorPoints[maxIndex]

        newHalfEdge=toDC.halfEdge(edge2.origin,None)
        newHalfEdge.destination=edge1.origin
        newHalfEdge.prev=edge2.prev
        newHalfEdge.next=edge1

        otherHalfEdge=toDC.halfEdge(edge1.origin,None)
        otherHalfEdge.destination=edge2.origin
        otherHalfEdge.prev=edge1.prev
        otherHalfEdge.next=edge2

        newHalfEdge.twin=otherHalfEdge
        otherHalfEdge.twin=newHalfEdge

        edge1.prev.next=otherHalfEdge
        edge1.prev=newHalfEdge

        edge2.prev.next=newHalfEdge
        edge2.prev=otherHalfEdge

        triangleList1=simpleTriangulation(edge1)
        triangleList2=simpleTriangulation(edge2)

    triangleList2.extend(triangleList1)
    return triangleList2

def isLeft(p1,p2,p3):
    """Given three points, we determine if p3 is on the left of the directed line given by points p1 and p2.
    
    Output: Signed distance between p3 and the line.
    When output>0, p3 is on the left of the line."""
    return ((p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x))/math.sqrt((p2.y-p1.y)**2+(p2.x-p1.x)**2)

def isInPolygon(point,pointList):
    """We use here a pointList to represent a polygon, with pointList[0] being pointList[-1].
    
    Return 0 only when point is outside. If the point is on the boundary, also return 1."""
    # We count the winding number.
    windingNumber=0

    for i in range(len(pointList)-1):
        if pointList[i].x==point.x and pointList[i].y==point.y: return 1
        if pointList[i+1].x==point.x and pointList[i+1].y==point.y: return 1
        if pointList[i].y==point.y==pointList[i+1].y and (point.x-pointList[i].x)*(point.x-pointList[i+1].x)<=0: return 1
        if pointList[i].y<=point.y:
            if pointList[i+1].y>point.y:
                if isLeft(pointList[i],pointList[i+1],point)>0: windingNumber+=1
                if isLeft(pointList[i],pointList[i+1],point)==0: return 1
        else:
            if pointList[i+1].y<=point.y:
                if isLeft(pointList[i],pointList[i+1],point)<0: windingNumber-=1
                if isLeft(pointList[i],pointList[i+1],point)==0: return 1
    
    return windingNumber