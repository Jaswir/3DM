def ChooseB (degree, x, y, z):
	if degree == 0:
		A = np.matrix([[1]])
		return A
	if degree == 1:
		A = np.matrix([[1],  [x], [y], [z], [x*y], [x*z], [y*z], [x*y*z]])
		return A
	if degree == 2:
		A = np.matrix([[1], [x], [y], [z], [x**2], [y**2], [z**2], [x*y], [x*z], [y*z], [x**2*y], [x**2*z], [x*y**2], [x*z**2], [y**2*z], [y*z**2], [x*y*z], [x**2*y*z], [x*y**2*z], [x*y*z**2], [x**2*y**2], [x**2*z**2], [y**2*z**2], [x**2*y**2*z], [x**2*y*z**2], [x*y**2*z**2], [x**2*y**2*z**2]])
		return A
print (ChooseB(2, 32, 52, 92))