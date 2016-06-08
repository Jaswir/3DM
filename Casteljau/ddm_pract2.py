# Place this file in your addons directory: e.g. "%appdata%\Roaming\Blender Foundation\Blender\2.77\scripts\addons", or install it as an addon using the Blender UI. Do not forget to enable the addon in the preferences panel.

bl_info = {
    "name": "DDM Practicum 2",
    "description": "Creates a subdivision surface according to Practicum 2.",
    "author": "Vanderfeesten",
    "version": (1, 0),
    "blender": (2, 77, 0),
    "location": "View3D > Add > Mesh",
    "warning": "", # used for warning icon and text in addons panel
    "wiki_url": "http://www.cs.uu.nl/docs/vakken/ddm/"
                "",
    "tracker_url": "",
    "support": "COMMUNITY",
    "category": "Add Mesh"
    }

import bpy
import random
import mathutils
from scipy.spatial import ConvexHull

from bpy.props import *

# DDM PRACTICUM 2
# ---------------

# See the assignment for details
# http://www.cs.uu.nl/docs/vakken/ddm/

# Check for outer vertex
def IsOutervertex(x, y, n):
    if y == 0 or y == (n - 1) or x == 0 or x == (n - 1):
        return True
    return False


#Visualize boundaries of 3D Grid
def makeMaterial(name, diffuse, specular, alpha):
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = diffuse
    mat.diffuse_shader = 'LAMBERT'
    mat.diffuse_intensity = 1.0
    return mat

def setMaterial(ob, mat):
    me = ob.data
    me.materials.append(mat)

# Order vertices first increasing in x, then in y
def ControlMesh(n, length): 
    vertices = []
    step = length/n
    for y in range (0, n):
        for x in range (0, n):
            z = 0
            if not(IsOutervertex(x, y, n)):
                z =  random.random()*step*2
            vertices.append(mathutils.Vector((y*step, x*step, z)))         
    return vertices


#Creates a mesh from vertices and triangles
#Source: https://wiki.blender.org/index.php/Dev:Py/Scripts/Cookbook/Code_snippets/Three_ways_to_create_objects  
def createMeshFromData(name, origin, vertices, faces):
    # Create mesh and object
    me = bpy.data.meshes.new(name+'Mesh')
    ob = bpy.data.objects.new(name, me)
    ob.location = origin
    ob.show_name = True
 
    # Link object to scene and make active
    scn = bpy.context.scene
    scn.objects.link(ob)
    scn.objects.active = ob
    ob.select = True
 
    # Create mesh from given vertices, faces.
    me.from_pydata(vertices, [], faces)
    # Update mesh with new data
    me.update()    

    return ob

#Given n, generates indices for faces.
#Where a face is a quad and has the format: [rightTopcorner, rightBotcorner, leftBotcorner, leftTopcorner]
def CreateFacesFromMesh(n):
    faces = []
    for y in range (0, n - 1):
        for x in range (0, n - 1):
            rightTopcorner = x + n + y * n + 1
            rightBotcorner = x + y * n + 1
            leftBotcorner  = x + y * n
            leftTopcorner  = x + n + y * n
            faces.append([rightTopcorner, rightBotcorner, leftBotcorner, leftTopcorner])
    return faces

#Displays the mesh on the screen
#also prints the vertices of the mesh
def ShowMesh(vertices, n, name):
    faces = CreateFacesFromMesh(n)
    createMeshFromData(name, (0,0,0), vertices, faces)
    #print(vertices)
    
#Given a list of vertices, amount of sqrt(vertices) and amount of subdivisions
#computes the appropriate DeCasteljau vertices. 
def DeCasteljau(A, n, s):

    #Returns A if the amount of subdivision is 0
    if s is 0:
        return A

    #Used to compute all the cas values
    stepSize = 0.5/(s+1)
    totalCasValues = int(1/stepSize)

    casCols = []
    for i in range(0,n):
        
        # Fills 1 column
        column = []
        for j in range(0, n):
             column.append(A[i+j*n])
 
        for k in range (0, totalCasValues+1):
            #Computes all the cas points for a  column  of the surface 
            casValue = k * stepSize

            #Holds the points for castelJau algorithm on a curve
            casPoints = []   
            result = column        
            length = n
                                        
            #Calculates the final point on the curve using casValue
            while(length != 1):
                casPoints = []
                for l in range(0, length-1):
                    casPoints.append(result[l]*(1-casValue) + result[l+1]*casValue)
                result = casPoints
                length -= 1
            casCols.append(result[0])

    

    #Using the result from the previous DeCasteljau for the vertical dimension,
    #DeCasteljau for the horizontal dimension is computed, resulting in 2D computed DeCastelJau.
    DeCasteljau2D = []
    for i in range(0, totalCasValues+1):
       
        #Fill 1 row
        row = []
        for j in range(0, n):
            row.append(casCols[i+j*(totalCasValues+1)])

   
        for k in range (0, totalCasValues+1):
            #Computes all the cas points for a row of the surface  
            casValue = k * stepSize

            #Holds the points for DeCasteljau algorithm on a curve
            casPoints = []   
            result = row
            length = n

            #Calculates the final point on the curve using casValue
            while(length != 1):
                casPoints = []
                for l in range(0, length-1):
                    casPoints.append(result[l]*(1-casValue) + result[l+1]*casValue)
                result = casPoints
                length -= 1
            DeCasteljau2D.append(result[0])

    return DeCasteljau2D
    
#Checks whether a point is inside a Rectangle (3_Dimensional)
#Type, indicates whether the face is BottomTop , Side or FrontBack
def pointInRect_3D(leftBottom, leftTop, rightBottom, rightTop, Point, typeFace):

    #If there's no intersection point, return False
    if(Point == None):
        return False

    i = None
    j = None

    #Checks which dimensions to choose 0: x, 1:y, 2: z
    if(typeFace == "Side"):
        i = 0 
        j = 2

    if(typeFace == "FrontBack"):
        i = 1
        j = 2


    if(typeFace == "BottomTop"):
        i = 1
        j = 0

    #Checks whether point is inside the rect using the proper dimensions
    if(Point[i] >= leftBottom[i] and Point[i] <= rightBottom[i] and Point[j] >= leftBottom[j] and Point[j] <= leftTop[j]):
        return True
    else:
        return False
    return True


#Returns boolean that indicates whether a line intersects a bounding box, given a list of vertices 
def LineIntersectsBoundingBox(A, p1, p2):
    #Gets boundary coördinates
    xList, yList, zList = [], [], []

    for i in range(0, len(A)):
        xList.append(A[i][0])
        yList.append(A[i][1])
        zList.append(A[i][2])

    xmin = min(xList)
    xmax = max(xList)
    ymin = min(yList)
    ymax = max(yList)
    zmin = min(zList)
    zmax = max(zList)

    #Gets 8 (axis-aligned bounding box) AABB points
    #When you take the rectangle on the left and right side you can bound all points.
    leftLeftTop = mathutils.Vector((xmin ,ymin, zmax)) 
    leftLeftBot = mathutils.Vector((xmin  ,ymin, zmin))
    leftRightTop = mathutils.Vector((xmax  ,ymin, zmax))
    leftRightBot = mathutils.Vector((xmax ,ymin, zmin))

    rightLeftBot = mathutils.Vector((xmin ,ymax, zmin))
    rightLeftTop = mathutils.Vector((xmin ,ymax, zmax)) 
    rightRightBot = mathutils.Vector((xmax ,ymax, zmin))
    rightRightTop = mathutils.Vector((xmax ,ymax, zmax))

    #Sets up preparations for plane-line intersection. 
    #3 planes of the AABB share the same point on plane
    #also opposing planes share the same normal
    #intersection algorithm
    normalOfsides = mathutils.geometry.normal(leftLeftTop, leftLeftBot, leftRightTop)
    normalOfBotandTop = mathutils.geometry.normal(leftLeftTop, rightRightTop, leftRightTop)
    normalOfFrontandBack = mathutils.geometry.normal(leftRightBot, leftRightTop, rightRightTop)
    pointOnBackLeftBotPlane = leftLeftBot
    pointOnFrontRightTopPlane = rightRightTop

    line_a = mathutils.Vector(p1)
    line_b = mathutils.Vector(p2)

    #Performs intersection for all six planes of the AABB
    intersections = []

    intersectionLeft = mathutils.geometry.intersect_line_plane(line_a, line_b, pointOnBackLeftBotPlane, normalOfsides)
    intersections.append(pointInRect_3D(leftLeftBot, leftLeftTop, leftRightBot, leftRightTop, intersectionLeft, "Side"))
    
    intersectionRight = mathutils.geometry.intersect_line_plane(line_a, line_b, pointOnFrontRightTopPlane, normalOfsides)
    intersections.append(pointInRect_3D(rightLeftBot, rightLeftTop, rightRightBot, rightRightTop, intersectionRight, "Side"))

    intersectionFront = mathutils.geometry.intersect_line_plane(line_a, line_b, pointOnFrontRightTopPlane, normalOfFrontandBack)
    intersections.append(pointInRect_3D(leftLeftBot, leftLeftTop, rightLeftBot, rightLeftTop, intersectionFront, "FrontBack"))

    intersectionBack = mathutils.geometry.intersect_line_plane(line_a, line_b, pointOnBackLeftBotPlane, normalOfFrontandBack)
    intersections.append(pointInRect_3D(leftRightBot, leftRightTop, rightRightBot, rightRightTop, intersectionBack, "FrontBack"))

    intersectionTop = mathutils.geometry.intersect_line_plane(line_a, line_b, pointOnFrontRightTopPlane, normalOfBotandTop)
    intersections.append(pointInRect_3D(leftLeftTop, leftRightTop, rightLeftTop, rightRightTop, intersectionTop, "BottomTop"))

    intersectionBot = mathutils.geometry.intersect_line_plane(line_a, line_b, pointOnBackLeftBotPlane, normalOfBotandTop)
    intersections.append(pointInRect_3D(leftLeftBot, leftRightBot, rightLeftBot, rightRightBot, intersectionBot, "BottomTop"))


    #Returns true if line intersects AABB, otherwise False
    for i in range(0,len(intersections)):
        if(intersections[i]):
            return True
    return False


def DeCasteljauSplit(A, n, p1, p2, e):

    beginPoints, endPoints = [], []
    casCols = []
    casValue = 0.5

    for i in range(0,n):
        
        # Fills 1 column
        column = []
        for j in range(0, n):
             column.append(A[j+i*n])
        
        #Holds the points for castelJau algorithm on a curve
        casPoints = []   
        result = column        
        length = n
        beginPoints, endPoints = [], []
                                            
        #Calculates the final point on the curve using casValue
        while(length != 1):
            casPoints = []
            for l in range(0, length-1):
                casPoints.append(result[l]*(1-casValue) + result[l+1]*casValue)
            beginPoints.append(result[0])
            if len(result) != 1:
                endPoints.append(result[len(result)-1])
            result = casPoints
            length -= 1
        beginPoints.append(result[0])

        points = beginPoints + list(reversed(endPoints))
        casCols = casCols + points


    #Using the result from the previous DeCasteljau for the vertical dimension,
    #DeCasteljau for the horizontal dimension is computed, resulting in 2D computed DeCastelJau.
    DeCasteljau2D = []
    newN = 2*n - 1
    for i in range(0, newN):
       
        #Fill 1 row
        row = []
        for j in range(0, n):
            row.append(casCols[i+j*(newN)])

        #Holds the points for DeCasteljau algorithm on a curve
        casPoints = []   
        result = row
        length = n
        beginPoints, endPoints = [], []

        #Calculates the final point on the curve using casValue
        while(length != 1):
            casPoints = []
            for l in range(0, length-1):
                casPoints.append(result[l]*(1-casValue) + result[l+1]*casValue)
            beginPoints.append(result[0])
            if len(result) != 1:
                endPoints.append(result[len(result)-1])
            result = casPoints
            length -= 1
        beginPoints.append(result[0])

        points = beginPoints + list(reversed(endPoints))
        DeCasteljau2D = DeCasteljau2D + points

    
   
    
    #Check whether shortest distance from middle point to line is smaller than e, if so returns true
    #See: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    indx = int((len(DeCasteljau2D)-1)/2)
    point = DeCasteljau2D[indx]
    p1Vec = mathutils.Vector(p1)
    p2Vec = mathutils.Vector(p2)
    shortestDistanceFromThisPointToLine = (((point - p1Vec).cross(point - p2Vec)).magnitude)/((p1Vec - p2Vec).magnitude)
    if(shortestDistanceFromThisPointToLine < e):
        return True
        

    size = 0.05
    for index, p in enumerate(DeCasteljau2D):
		#LB
	    red = makeMaterial('Red',(1,0,0),(1,1,1),1)
	    origin = p
	    bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
	    bpy.ops.transform.resize(value=(size, size, size))
	    setMaterial(bpy.context.object, red)
    
    #Form 4 surfaces, assign all surfaces their portion of the control points
    #Go in recursion on these
    splittedLength = int((newN-1)/2+1)
    leftBotCornerSurfacePoints = []
    rightBotCornerSurfacePoints = []
    leftTopCornerSurfacePoints  = []
    rightTopCornerSurfacePoints = []
    for y in range (0, splittedLength):
        for x in range (0, splittedLength):
            leftBotCornerSurfacePoints.append(DeCasteljau2D[x+ y*newN])
            rightBotCornerSurfacePoints.append(DeCasteljau2D[splittedLength -1 + x + y*newN])
            leftTopCornerSurfacePoints.append(DeCasteljau2D[newN*2 + x + (y+(n-3))*newN])
            rightTopCornerSurfacePoints.append(DeCasteljau2D[newN*2 + splittedLength -1 + x + (y+(n-3))*newN])
            
 
    '''
    LineIntersect(leftBotCornerSurfacePoints, splittedLength, p1, p2, e)
    LineIntersect(rightBotCornerSurfacePoints, splittedLength, p1, p2, e)
    LineIntersect(leftTopCornerSurfacePoints,  splittedLength , p1, p2, e)
    LineIntersect(rightTopCornerSurfacePoints,  splittedLength, p1, p2, e)
    ''' 
   	
    size = 0.2
    for index in range(splittedLength*splittedLength):
    	#LB
        green = makeMaterial('Green',(0,1,0),(1,1,1),1)
        origin = leftBotCornerSurfacePoints[index]
        bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
        bpy.ops.transform.resize(value=(size, size, size))
        setMaterial(bpy.context.object, green)
        
        #RB
        blue = makeMaterial('blue',(0,0,1),(1,1,1),1)
        origin = rightBotCornerSurfacePoints[index]
        bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
        bpy.ops.transform.resize(value=(size, size, size))
        setMaterial(bpy.context.object, blue)

        
        #LT
        purple = makeMaterial('purple',(0.7,0,1),(1,1,1),1)
        origin = leftTopCornerSurfacePoints[index]
        bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
        bpy.ops.transform.resize(value=(size, size, size))
        setMaterial(bpy.context.object, purple)

        #RT
        red = makeMaterial('red',(1,0,0),(1,1,1),1)
        origin = rightTopCornerSurfacePoints[index]
        bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
        bpy.ops.transform.resize(value=(size, size, size))
        setMaterial(bpy.context.object, red)
    	
     
	
    
    
#Assume inputted array is filled with Vectors
def LineIntersect(A, n, p1, p2, e):

    #Compute bounding box of surface
    if not(LineIntersectsBoundingBox(A, p1, p2)):
        return False
    else:
        return DeCasteljauSplit(A, n, p1, p2, e)
            
        


          


    
def main(operator, context):
    
    n = 5
    length = 10
    s = 10
    #DeCasteljau needs different n for faces, which can be computed as follows.
    nFaceCas = int(1+ 1/(0.5/(s+1)))
    
    A = ControlMesh(n, length)
    B = DeCasteljau(A, n, s)

    ShowMesh(A, n, "display")
    #ShowMesh(B, nFaceCas, "CastelJau")
    p1 = (3,4,3)
    p2 = (3,4,7)

    print(LineIntersect(A, n, p1, p2, 0.01))

# BLENDER UI
# ----------

class MESH_OT_primitive_ddm_surface_add(bpy.types.Operator):
    '''Add a subdivision surface'''
    bl_idname = "mesh.primitive_ddm_surface_add"
    bl_label = "Add subdivision surface"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        main(self, context)
        return {'FINISHED'}

def menu_func(self, context):
    self.layout.operator("mesh.primitive_ddm_surface_add",
        text="DDM Subdivision Surface",
        icon='MESH_PLANE')

def register():
   bpy.utils.register_module(__name__)
   bpy.types.INFO_MT_mesh_add.prepend(menu_func)

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_mesh_add.remove(menu_func)

if __name__ == "__main__":
    register()