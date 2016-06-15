# Either run the script from the internal Blender script editor or execute the script on the console. If you are using an external editor, use the following commands on the blender console to run the script (this always takes the new version)
#
# filename = "/DDM/ddm_pract3.py"
# exec(compile(open(filename).read(), filename, 'exec'), globals(), locals())

import bpy
import math
import bmesh
import numpy as np
import scipy
from mathutils import *


# Runs the precalculation process, storing the precomputation data with the object that is being deformed
def Precompute(source_object):




	# Get a BMesh representation
	bm = bmesh.new()              # create an empty BMesh
	bm.from_mesh(source_object.data)   # fill it in from a Mesh


	#Turn A into A' by removing constraint(CONST) columns.
	#Instead of removing constraint columns, we omit adding these columns to A altogether.

	#Get const column indices
	CONST = get_handles(source_object)[0][0]


	#The A' matrix 
	APrime = []
	idx = 0
	CONSTidx = 0;
	for vertex in bm.verts:


		validColumn = idx != CONST[CONSTidx]
		
		#Check whether current vertex idx is a CONST index.
		#Every time this is false we add the column to APrime
		if(idx != ):
			#Holds weights for 1 one-ring
			column = []
			#Get oneRing neighbours of vertex
			one_ringNeighbours = [edge.other_vert(vertex) for edge in vertex.link_edges]


			for neighbour in one_ringNeighbours:
				#Get oneRing neighbours of neighbours
				neighboursNeighbours = [edge.other_vert(neighbour) for edge in neighbour.link_edges]
				#Find matching neighbours
				matchingNeighbours =  list(set(one_ringNeighbours) & set(neighboursNeighbours))

				#Use these to compute weights, see https://in.answers.yahoo.com/question/index?qid=20110530115210AA2eAW1
				aAlpha = neighbour.co - matchingNeighbours[0].co
				aBeta  = neighbour.co - matchingNeighbours[1].co
				bAlpha = vertex.co - matchingNeighbours[0].co
				bBeta  = vertex.co - matchingNeighbours[1].co

				'''MIGHT NEED TO BE ADAPTED IN THE FUTURE IF THINGS LOOK WRONG'''
				tanAlpha = ((aAlpha.cross(bAlpha)).magnitude) / (aAlpha.dot(bAlpha))
				tanBeta = ((aBeta.cross(bBeta)).magnitude) / (aBeta.dot(bBeta))
				cotAlpha = 1/tanAlpha
				cotBeta = 1/tanBeta

				Wiv = 0.5*(cotAlpha + cotBeta)
				#AMIR: Some people asked what to do with negative cotangent weight, 
				#because they commonly take sqrt(w_ij) to put into the ||Ax-b||^2 expression. 
				#The weights should not overly negative in reasonable triangles, 
				#so try and use w_{ij}={small positive epsilon} as a cheap workaround. For instance w_{ij}=10e-3.
				if (Wiv < 0):
					Wiv = 10e-3

				column.append(Wiv)

			column = np.array(column)
			column = np.sqrt(column)

			APrime.append(column)

		#Otherwise we update the CONSTidx,
		#so it will grab the next index in the following iteration
		#when the last index of CONSTidx has been reached
		#all constraint columns are taken account of
		else:
			CONSTidx ++

		#Prepare idx for the iteration
		idx ++

	APrime = np.transpose(np.array(APrime))









	
	#source_object.data['precomputed_data']  = Matrix()


# Runs As Rigid As Possible deformation on the mesh M, using a list of handles given by H. A handle is a list of vertices tupled with a transform matrix which might be rigid (identity)
#def ARAP(source_mesh, deformed_mesh, H):
	

def main():

	'''
	A = []
	a = np.array([1,2,3])
	b = np.array([4,5,6])
	c = np.array([7,8,9])

	A.append(a)
	A.append(b)
	A.append(c)

	A = np.array(A)
	print(A)
	print("\n")


	A = scipy.delete(A , 1  ,1)
	print(A)
	'''

	# TODO: Check for an existing deformed mesh, if so use that as an iteration, if not use a mesh named 'source' as the initial mesh.
	source = bpy.data.objects['source']
	
	# TODO: Precompute A'^T * A if the data is dirty or does not exist, and store it with the source object
	Precompute(source)
	
	# TODO: Perform As Rigid As Possible deformation on the source object in the first iteration, and on a deformed object if it exists
	#ARAP(source.data, get_deformed_object(source).data, get_handles(source))


	
	
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
		
	return T * Minv

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
			
			# Interpret the transform as a bounding box for selecting the handles
			handle_vertices = get_handle_vertices(source.data.vertices, bounding_box_transform, mesh_transform)
			
			# If a destination box exists
			handle_dest_name = handle_name + '_dest'
			if bpy.data.objects.get(handle_dest_name) is not None:
				
				bounding_box_dest_transform = get_transform_of_object(handle_dest_name)
				
				result.append( (handle_vertices, get_relative_transform(bounding_box_transform, bounding_box_dest_transform) ) ) 
				
			else:
			
				# It is a rigid handle
				m = Matrix()
				m.identity()
				result.append( (handle_vertices, m) )
			
	return result
	
main()