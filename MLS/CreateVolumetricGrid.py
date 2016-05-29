import bpy
import mathutils
import mcubes
import numpy as np


def main(context):

#1. Dimension of bounding box of object (USED TO SAMPLE 3D GRID)  #
##################################################    #


    #1. HERE THE BOUNDARIES OF THE OBJECT ARE COMPUTED ---
    name = context.active_object.data.name


    xDimensionsHalve = bpy.data.objects[name].dimensions.x/2  
    xMin = 0 - xDimensionsHalve
    xMax = 0 + xDimensionsHalve

    yDimensionsHalve = bpy.data.objects[name].dimensions.y/2
    yMin = 0 - yDimensionsHalve
    yMax = 0 + yDimensionsHalve

    zDimensionsHalve = bpy.data.objects[name].dimensions.z/2
    zMin = 0 - zDimensionsHalve
    zMax = 0 + zDimensionsHalve



#2. COMPUTE MLS FUNCTION ##########################   #
##################################################    #

    #0. Define Help functions:
    #
    #Wendland weight function
    def Wendland ( r , h ):
         return (1 - r/h)**4*(4*r/h+1)

    #b^T for degree 0, 1 or 2.
    def ChooseB (degree, point):
        x = point[0]
        y = point[1]
        z = point[2]
        if degree == 0:
             A = np.matrix([[1]])        
        elif degree == 1:
             A = np.matrix([[1],  [x], [y], [z], [x*y], [x*z], [y*z], [x*y*z]])
        elif degree == 2:
             A = np.matrix([[1], [x], [y], [z], [x**2], [y**2], [z**2], [x*y], [x*z], [y*z], [x**2*y], [x**2*z], [x*y**2], [x*z**2], [y**2*z], [y*z**2], [x*y*z], [x**2*y*z], [x*y**2*z], [x*y*z**2], [x**2*y**2], [x**2*z**2], [y**2*z**2], [x**2*y**2*z], [x**2*y*z**2], [x*y**2*z**2], [x**2*y**2*z**2]])
        return np.transpose(A)

    #Checks whether a point P is within a range of the given point C
    #which in 3D means is within the sphere with C as its center and the specified
    #, range as its radius
    def insideSphere(C, R, P):
        return ( P[0]- C[0] ) ** 2 + (P[1]- C[1]) ** 2 + (P[2]-C[2]) ** 2 < R**2 

    #Gets the appropriate length that matches to a degree k
    def matchingLengthofDegreeK(k):
        if k == 0:
             l = 1
        elif k == 1:
             l = 8
        elif k == 2:
             l = 27
        return l  

    #Can be used for visualization purposes
    def makeMaterial(name, diffuse, specular, alpha):
        mat = bpy.data.materials.new(name)
        mat.diffuse_color = diffuse
        mat.diffuse_shader = 'LAMBERT'
        mat.diffuse_intensity = 1.0
        return mat

    def setMaterial(ob, mat):
        me = ob.data
        me.materials.append(mat)

    #1. SETUP EPSILON, POINTCLOUD POINTS (FIRST n POINTS), POINTCLOUD POINT NORMALS and N (NECESSARY VARIABLES FOR COMPUTING CONSTRAINT POINTS), POINTCLOUD CENTER --- 
    #
    #Computes N (amount of points in pointCloud)
    PointCloud_points = []
    PointCloud_points  = context.active_object.data.vertices
    N = len(PointCloud_points )

    #Computes Normals of point cloud points
    PointCloud_point_normals = []
    for n in context.active_object['vertex_normal_list']:
        PointCloud_point_normals.append(np.array([n[0], n[1], n[2]]))


    #Fixes an ε value, for instance ε = 0.01 times the diagonal of the bounding box to the object. EXPERIMENTABLE PARAMETER
    v1 =  np.array((xMin, yMin,zMax))
    v2 =  np.array((xMax, yMax,zMin))
    diagonalOfBoundingBox = np.linalg.norm(v1-v2)
    epsilon = 0.01 * diagonalOfBoundingBox

    #2. IMPLEMENT SPATIAL INDEX: KD-TREE ---
    #
    #Computes a Spatial Index: KD-TREE from the PointCloud_points
    #(for faster nearest neighborhood calculations)
    kd = mathutils.kdtree.KDTree(N)
    for index, point in enumerate(PointCloud_points):
        kd.insert(point.co, index)
        
    #Must have been called before using any of the kd.find methods.
    kd.balance()


    #Compute the origin of the axis aligned bounding box
    origin = mathutils.Vector((0,0,0))
    for index, point in enumerate(PointCloud_points):
         origin += point.co

    origin /= N 
    polyDegree = 0
    def MLS(x , y ,z):
                 P = np.array([x ,  y,  z])
                 #Get closest points Pi within wendland radius from P. EXPERIMENTABLE PARAMETER
                 wendlandRadius = diagonalOfBoundingBox/20
                 pointsWithinRange = kd.find_range(P, wendlandRadius)
                 numpyPointsWithinRange = []
                 redPoints = []
                 redPointsEpsilonValues = []
                 greenPoints = []
                 greenPointsEpsilonValues = []
                 for (co,index, dist) in pointsWithinRange:
                     Ni = PointCloud_point_normals[index]
                     pos = co
                     PiNumpy = np.array((pos[0], pos[1], pos[2]))
                     epsilonBackup = epsilon
                     result = PiNumpy + epsilonBackup*Ni  
                     #Range to look for closest points from Pi+N. EXPERIMENTABLE PARAMETER              
                     Range = diagonalOfBoundingBox/20
                     Pi_Not_Closest = True
                     #Change datatype of Pi, to make computations easier
                     Pi = pos 
                     #Check whether Pi is the closest point to Pi+N
                     while Pi_Not_Closest :
                         #Calculates distances between Pi+N = result and the points that are within a range distance from Pi+N
                         dist, closestPos = min([(np.linalg.norm(result - co), co) for (co,index, dist) in kd.find_range(result, Range)])
  
                         #If Pi is not the closest point to Pi+N, 
                         # divide ε by 2 and recompute pi+N until this is the case
                         if closestPos != Pi:
                             epsilonBackup /= 2
                             result = Pi + epsilonBackup * Ni
                         else:
                             Pi_Not_Closest = False
                             break    

                     #check whether the green and red point are within wendlandRadius of P  
                     if insideSphere(C = P, R = wendlandRadius, P = result - 2*(epsilonBackup * Ni)):
                         greenPoints.append(result - 2*(epsilonBackup * Ni))
                         greenPointsEpsilonValues.append(- epsilonBackup)

                     if insideSphere(C = P, R = wendlandRadius ,P = result):                      
                         redPoints.append((result))
                         redPointsEpsilonValues.append(epsilonBackup)
                     numpyPointsWithinRange.append(PiNumpy)

                 #Create lists of the Pi's and Di's belonging to P
                 PointsPi = numpyPointsWithinRange + redPoints + greenPoints
                 PointsDi = [0]*len(pointsWithinRange) + redPointsEpsilonValues + greenPointsEpsilonValues

                 #sqrt(ti) = sqrt(theta(|P-Pi|)) (see MLS reference sheet: http://www.cs.uu.nl/docs/vakken/ddm/MLS%20Reference%20Sheet.pdf)
                 sqrtTiValues = [] 
                 #smoothing value for wendland weight function. EXPERIMENTABLE PARAMETER
                 h = 1
                 for Pi in PointsPi: 
                     value = np.sqrt(Wendland(np.linalg.norm(P-Pi), h))        
                     sqrtTiValues.append(value)
                 # A is a one-column matrix filled by sqrtTi multiplied with b^T
                 A = np.empty((0,matchingLengthofDegreeK(polyDegree)))
                 for index, sqrtTi in enumerate(sqrtTiValues):
                     #Here we choose B for a degree 0, 1 or 2. EXPERIMENTABLE PARAMETER
                     bT = ChooseB(polyDegree, P)
                     value = sqrtTi*bT
                     A = np.insert(A, index, value, 0)

                    
                 # r = sqrtTi_x_Di
                 r = np.empty((0, len(sqrtTiValues)))
                 for index, sqrtTi in enumerate(sqrtTiValues):
                     r = np.insert(r, index, sqrtTi * PointsDi[index])
                 

                 if(len(A) != 0):
                     A_T = np.transpose(A)
                     A_T_x_A = np.dot(A_T, A)        
                     A_T_x_r = np.dot(A_T, r)
                     #a = (A^T*A)^-1 * A^T*r
                     a = np.dot(np.linalg.inv(A_T_x_A), A_T_x_r)
                     #Finally add a to the samples    
                     return np.dot(ChooseB(polyDegree, [x,y,z]), a)

                 else:
                     return 10000

    def f(x, y, z):
        return MLS(x, y, z)


    #4. INPUTS SAMPLED GRID TO MARCHING CUBES PLUGIN ---
    #
    lowerLeft = origin - mathutils.Vector((diagonalOfBoundingBox/2*0.9, diagonalOfBoundingBox/2*0.9, diagonalOfBoundingBox/2*0.9))
    upperRight = origin + mathutils.Vector((diagonalOfBoundingBox/2*1.1, diagonalOfBoundingBox/2*1.1, diagonalOfBoundingBox/2*1.1))

    vertices, triangles = mcubes.marching_cubes_func((lowerLeft[0], lowerLeft[1], lowerLeft[2]),(upperRight[0], upperRight[1], upperRight[2]), 100, 100, 100, f, 0)

    # Export the result
    mcubes.export_mesh(vertices, triangles, "C:\\Users\\jaswir\\Documents\\GameTechnology\\3DM\\3DM_Practical1\\DAE_Files\\Bunnyk0.dae", "Bunny_k0")


class CreateVolumetricGrid(bpy.types.Operator):
    """Tooltip"""
    #labelname Create Volume, because it sounds cool. 
    bl_idname = "object.simple_operator"
    bl_label = "Create Volume"

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        main(context)
        return {'FINISHED'}

#Register the class to blender
def register():
    bpy.utils.register_class(CreateVolumetricGrid)


def unregister():
    bpy.utils.unregister_class(CreateVolumetricGrid)


if __name__ == "__main__":
    register()
