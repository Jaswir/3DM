# Either run the script from the internal Blender script editor or execute the script on the console. If you are using an external editor, use the following commands on the blender console to run the script (this always takes the new version)
#
# filename = "/DDM/ddm_pract3.py"
# exec(compile(open(filename).read(), filename, 'exec'), globals(), locals())

import bpy
import bmesh
import numpy as np




# Runs the precalculation process, storing the precomputation data with the object that is being deformed
def Precompute(source_object):

	# Get a BMesh representation
	bm = bmesh.new()              # create an empty BMesh
	bm.from_mesh(source_object.data)   # fill it in from a Mesh

	#Questions: 
	#-1 Deformed and original mesh??
	#-2 SVD - fullmatrices option

	'''NOTES
	If you use the indices, you probably don't have to take account of the changing positions as you can call the smae indices.
	'''
	
	#Computes Local Step 

	for vertex in bm.verts:
		#Composes Matrices P and Q 
		one_ringNeighbours = []
		for edge in vertex.link_edges:
			neighbour = edge.other_vert(vertex)
			one_ringNeighbours.append(vertex.co - neighbour.co)
		one_ringNeighbours = np.array(one_ringNeighbours)

		
		pT = np.transpose(one_ringNeighbours)
		
		Q = ? 

		#Compose correlation matrix S3×3 = P^T Q.
		S3x3 = np.dot(pT, Q)  

		#Decompose S = UΣvT using singular value decomposition, 
		#and composes the closest rigid transformation Rv = UvT
		U, s, v = np.linalg.svd(S3x3)
		vT = np.transpose(v)
		Rv = np.dot(U, vT)

		detRv = np.linalg.det(Rv)
		if(detRv == -1):
			




# Runs As Rigid As Possible deformation on the mesh M, using a list of handles given by H. A handle is a list of vertices tupled with a transform matrix which might be rigid (identity)
def ARAP(source_mesh, deformed_mesh, H):

	# TODO: remove this code and implement
	print(source_mesh.vertices)
	print(deformed_mesh.vertices)
	print(H)
	return source_mesh

def main():
	# TODO: Check for an existing deformed mesh, if so use that as an iteration, if not use a mesh named 'source' as the initial mesh.
	source = bpy.context.active_object
	
	# TODO: Precompute A'^T * A if the data is dirty or does not exist, and store it with the source object
	Precompute(source)
	
	# TODO: Perform As Rigid As Possible deformation on the source object in the first iteration, and on a deformed object if it exists
	#ARAP(source.data, get_deformed_object(source).data, get_handles(source))
	#get_handles(source)    
	
	
# BLENDER
# -------
	
def get_transform_of_object(name):
	return bpy.data.objects[name].matrix_basis

def get_mesh_vertices(name):
	return bpy.data.objects[name].data.vertices
	
# Finds the relative transform from matrix M to T
def get_relative_transform(M, T):
	Minv = M.copy()
	Minv.invert()
	return Minv * T

# Returns an object that can be used to store the deformed mesh
def get_deformed_object(source):
	
	name = 'deformed'
	
	# Create an object if it doesn't yet exist
	if bpy.data.objects.get(name) is None:
	
		# Create new mesh
		mesh = bpy.data.meshes.new(name)
	 
		# Create new object associated with the mesh
		ob_new = bpy.data.objects.new(name, mesh)
	 
		scn = bpy.context.scene
		scn.objects.link(ob_new)
		scn.objects.active = ob_new
	 
		# Copy data block from the old object into the new object
		ob_new.data = source.data.copy()
		ob_new.scale = source.scale
		ob_new.location = source.location
	
	return bpy.data.objects[name]
	
# Find the vertices within the bounding box by transforming them into the bounding box's local space and then checking on axis aligned bounds.
def get_handle_vertices(vertices, bounding_box_transform, mesh_transform):

	result = []

	# Fibd the transform into the bounding box's local space
	bounding_box_transform_inv = bounding_box_transform.copy()
	bounding_box_transform_inv.invert()
	
	# For each vertex, transform it to world space then to the bounding box local space and check if it is within the canonical cube x,y,z = [-1, 1]
	for i in range(len(vertices)):
		vprime = vertices[i].co.copy()
		vprime.resize_4d()
		vprime = bounding_box_transform_inv * mesh_transform * vprime
		
		x = vprime[0]
		y = vprime[1]
		z = vprime[2]
		
		if (-1 <= x) and (x <= 1) and (-1 <= y) and (y <= 1) and (-1 <= z) and (z <= 1):
			result.append(i)

	return result

# Returns a list of handles and their transforms
def get_handles(source):
	
	result = []
	
	mesh_transform = get_transform_of_object(source.name)
	
	# Only search up to (and not including) this number of handles
	max_handles = 10
	
	# For all numbered handles
	for i in range(max_handles):
	
		# Construct the handles representative name
		handle_name = 'handle_' + str(i)
		
		# If such a handle exists
		if bpy.data.objects.get(handle_name) is not None:
			
			# Find the extends of the aligned bounding box
			bounding_box_transform = get_transform_of_object(handle_name)
			print(bounding_box_transform)
			# Interpret the transform as a bounding box for selecting the handles
			handle_vertices = get_handle_vertices(source.data.vertices, bounding_box_transform, mesh_transform)
			
			# If a destination box exists
			handle_dest_name = handle_name + '_dest'
			if bpy.data.objects.get(handle_dest_name) is not None:
				
				bounding_box_dest_transform = get_transform_of_object(handle_dest_name)
				
				result.append( (handle_vertices, get_relative_transform(bounding_box_transform, bounding_box_dest_transform) ) ) 
				
			else:
			
				# It is a rigid handle
				result.append( (handle_vertices, Matrix([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]) ) )
			
	return result
	
main()