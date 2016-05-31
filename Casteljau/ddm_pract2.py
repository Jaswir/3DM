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
                z = random.uniform(step, length/2)
            size = 0.01
            vertices.append((y*step, x*step, z))
            red = makeMaterial('Red',(1,0,0),(1,1,1),1)
            origin = (y*step, x*step, z)
            bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
            bpy.ops.transform.resize(value=(size, size, size))
            setMaterial(bpy.context.object, red)
    return vertices

def ShowMesh(A, n):
    pass
    
def DeCasteljau(A, n, s):
    return []
    
def LineIntersect(A, n, p1, p2, e):
    return False
    
def main(operator, context):
    
    n = 10
    length = 1
    s = 3
    
    A = ControlMesh(n, length)
    B = DeCasteljau(A, n, s)
    ShowMesh(B, n)
    
    p1 = (1,2,3)
    p2 = (3,4,5)

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