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
    print(vertices)
    
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
    
def LineIntersect(A, n, p1, p2, e):
    return False
    
def main(operator, context):
    
    n = 5
    length = 10
    s = 10
    #DeCasteljau needs different n for faces, which can be computed as follows.
    nFaceCas = int(1+ 1/(0.5/(s+1)))
    
    A = ControlMesh(n, length)
    B = DeCasteljau(A, n, s)

    ShowMesh(A, n, "display")
    ShowMesh(B, nFaceCas, "CastelJau")
    
    p1 = (1,2,3)
    p2 = (3,4,5)

    #print(LineIntersect(A, n, p1, p2, 0.01))

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