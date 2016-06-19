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
import mathutils
from scipy.sparse import csr_matrix
import time


# Runs the precalculation process, storing the precomputation data with the object that is being deformed
def Precompute(source_object):
	
	# Get a BMesh representation
	bm = bmesh.new()              # create an empty BMesh
	bm.from_mesh(source_object.data)   # fill it in from a Mesh



	#Get const columns' indices
	CONST = get_handles(source_object)[0][0]
	
	amountOfcolumns = len(bm.verts)

	#Use KD-Tree to link vertices to indices in bm.verts
	size = len(bm.verts)
	verticesKD_Tree = mathutils.kdtree.KDTree(size)
	for i, v in enumerate(bm.verts):
		verticesKD_Tree.insert(v.co, i)

	verticesKD_Tree.balance()

	#Make new list, because bm.verts was being a bitch.
	vertices  = [vert.co for vert in bm.verts]
	
	#Finds the total amount of rows
	amountOfRows = 0
	for index, vertex in enumerate(bm.verts):
			one_ringNeighbours = [edge.other_vert(vertex) for edge in vertex.link_edges]
			amountOfRows += 1 + len(one_ringNeighbours)

	#Stores A
	A = np.zeros((amountOfRows, amountOfcolumns))

	#Used to count which row we are currently working at
	rowCount = 0

	for index, vertex in enumerate(bm.verts):

		#Get oneRing neighbours of vertex
		one_ringNeighbours = [edge.other_vert(vertex) for edge in vertex.link_edges]

		for neighbour in one_ringNeighbours:

			#Get oneRing neighbours of neighbour
			neighboursNeighbours = [edge.other_vert(neighbour) for edge in neighbour.link_edges]
			#Find matching neighbours
			matchingNeighbours =  list(set(one_ringNeighbours) & set(neighboursNeighbours))
	
			#Use these to compute weights, see https://in.answers.yahoo.com/question/index?qid=20110530115210AA2eAW1
			aAlpha = neighbour.co - matchingNeighbours[0].co
			aBeta  = neighbour.co - matchingNeighbours[1].co
			bAlpha = vertex.co - matchingNeighbours[0].co
			bBeta  = vertex.co - matchingNeighbours[1].co

			tanAlpha = ((aAlpha.cross(bAlpha)).magnitude) / (aAlpha.dot(bAlpha))
			tanBeta  = ((aBeta.cross(bBeta)).magnitude) / (aBeta.dot(bBeta))
			cotAlpha = 1/tanAlpha
			cotBeta  = 1/tanBeta 

			Wiv = 0.5*(cotAlpha + cotBeta)
			#AMIR: Some people asked what to do with negative cotangent weight, 
			#because they commonly take sqrt(w_ij) to put into the ||Ax-b||^2 expression. 
			#The weights should not overly negative in reasonable triangles, 
			#so try and use w_{ij}={small positive epsilon} as a cheap workaround. For instance w_{ij}=10e-3.
			if (Wiv < 0):
				Wiv = 10e-3

			Wiv = math.sqrt(Wiv)

			#Fills A
			vertexIdx    = verticesKD_Tree.find(vertex.co)[1]
			neighbourIdx = verticesKD_Tree.find(neighbour.co)[1]
			A[rowCount][vertexIdx]    = Wiv
			A[rowCount][neighbourIdx] = -Wiv

			#update row
			rowCount += 1
			
	#Stores oldA, for computing b'		
	oldA = A
	#Acquires A', by removing constraint columns from A.
	for columnIdx in reversed(CONST):
		A = np.delete(A,columnIdx,1)
	APrime = A

	#Computes A'TA' and prefactors the matrix
	APrime = np.matrix(APrime)
	APrimeTAPrime = APrime.T * APrime
	L = np.linalg.cholesky(APrimeTAPrime)

	# Saves things which will be used in further steps of the computation in dictionary
	source_object.data['precomputed_data'] = np.array(L)
	source_object.data['amountOfRows'] = amountOfRows
	source_object.data['A'] = np.array(A)
	
# Runs As Rigid As Possible deformation on the mesh M, using a list of handles given by H. A handle is a list of vertices tupled with a transform matrix which might be rigid (identity)
def ARAP(source_mesh, deformed_mesh, H):
	 # Get a BMesh representation
	bm = bmesh.new()              # create an empty BMesh
	bm.from_mesh(source_mesh)   # fill it in from a Mesh

	'''Handle existing deformed mesh'''

	#Computes Local Step 
	Rvs = []
	for vertex in bm.verts:
		#Composes Matrices P and Q 
		one_ringNeighbours = []
		for edge in vertex.link_edges:
			neighbour = edge.other_vert(vertex)
			one_ringNeighbours.append(vertex.co - neighbour.co)
		one_ringNeighbours = np.array(one_ringNeighbours)

		
		p = one_ringNeighbours
		pT = np.transpose(p)
		Q = p

		#Compose correlation matrix S3×3 = P^T Q.
		S3x3 = np.dot(pT, Q)  


		#Decompose S = UΣvT using singular value decomposition, 
		#and composes the closest rigid transformation Rv = UvT
		U, sigma, v = np.linalg.svd(S3x3)
		vT = np.transpose(v)
		Rv = np.dot(U, vT)


		#If det(Rv) = −1 (determinant), leading to reflection, instead compute Rv = UΣvT
		# where Σ'is an identity matrix, save for Σ'ii = -1 , where i is the index of 
		#the smallest diagonal (singular) value in the original Σ. . 
		#(flipping sign). For instance, if i = 3, you should use Σ' = diag[1, 1, −1].
		detRv = np.linalg.det(Rv)
		if(detRv == -1):
			#Find the index of the smallest diagonal vlaue in the original Σ.
			i =  np.argmin(sigma)
			#Computes Σ'
			s = [1,1,1]
			s[i] = -1
			sigmaTag = np.array([[s[0],0,0], [0,s[1],0], [0,0,s[2]]])
			#Computes reflection resistant  Rv
			Rv = np.dot(np.dot(U, sigmaTag), vT)
		print(Rv)

		Rvs.append(Rv)


		#Sets up b
		amountOfcolumns = 3
		amountOfRows = int(source_object.data['amountOfRows'])

		#Stores b
		b = np.zeros((amountOfRows, amountOfcolumns))

		#Use KD-Tree to link vertices to indices in bm.verts
		size = len(bm.verts)
		verticesKD_Tree = mathutils.kdtree.KDTree(size)
		for i, v in enumerate(bm.verts):
			verticesKD_Tree.insert(v.co, i)
		verticesKD_Tree.balance()

		#Make new list, because bm.verts was being a bitch.
		vertices  = [vert.co for vert in bm.verts]

		#Used to count which row we are currently working at
		rowCount = 0

		#Fills b
		for index, vertex in enumerate(bm.verts):

			#Get oneRing neighbours of vertex
			one_ringNeighbours = [edge.other_vert(vertex) for edge in vertex.link_edges]

			for neighbour in one_ringNeighbours:

				#Get oneRing neighbours of neighbour
				neighboursNeighbours = [edge.other_vert(neighbour) for edge in neighbour.link_edges]
				#Find matching neighbours
				matchingNeighbours =  list(set(one_ringNeighbours) & set(neighboursNeighbours))
		
				#Use these to compute weights, see https://in.answers.yahoo.com/question/index?qid=20110530115210AA2eAW1
				aAlpha = neighbour.co - matchingNeighbours[0].co
				aBeta  = neighbour.co - matchingNeighbours[1].co
				bAlpha = vertex.co - matchingNeighbours[0].co
				bBeta  = vertex.co - matchingNeighbours[1].co

				tanAlpha = ((aAlpha.cross(bAlpha)).magnitude) / (aAlpha.dot(bAlpha))
				tanBeta  = ((aBeta.cross(bBeta)).magnitude) / (aBeta.dot(bBeta))
				cotAlpha = 1/tanAlpha
				cotBeta  = 1/tanBeta 

				Wiv = 0.5*(cotAlpha + cotBeta)
				#AMIR: Some people asked what to do with negative cotangent weight, 
				#because they commonly take sqrt(w_ij) to put into the ||Ax-b||^2 expression. 
				#The weights should not overly negative in reasonable triangles, 
				#so try and use w_{ij}={small positive epsilon} as a cheap workaround. For instance w_{ij}=10e-3.
				if (Wiv < 0):
					Wiv = 10e-3

				#RowInput => sqrt(Wiv)(xv-xi)Rv
				vertexIdx = verticesKD_Tree.find(vertex.co)[1]
				rowInput  = math.sqrt(Wiv)*(vertex.co - neighbour.co)*Rvs[vertexIdx]
				
				b[rowCount]=  rowInput

				#update row
				rowCount += 1


		print("THUNDERSTRUCK YEAHYEAHYEAH")
		#b'= b - A (   0   )
		#          (X'Const)
		A = np.matrix(source_object.data['A'])
		#Zero..XPrimeConst = 
		#bPrime = b - A * Zero..XPrimeConst

   
	

def main():

	start = time.time()
	# TODO: Check for an existing deformed mesh, if so use that as an iteration, if not use a mesh named 'source' as the initial mesh.
	source = bpy.data.objects['source']
	
	# TODO: Precompute A'^T * A if the data is dirty or does not exist, and store it with the source object
	Precompute(source)
	
	# TODO: Perform As Rigid As Possible deformation on the source object in the first iteration, and on a deformed object if it exists
	ARAP(source.data, get_deformed_object(source).data, get_handles(source))

	print('It took {0:0.1f} seconds'.format(time.time() - start))


	
	
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