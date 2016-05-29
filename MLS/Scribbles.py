''' --- Commented Out ---
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


    #Visualizing axis aligned bounding box
    size = 0.05
    for i in range(8):
        x = None;
        y = None;
        z = None;
        if i >3:
            z = zMax
        else:
            z = zMin
        if i % 2 == 0:
            y = yMin
        else:
            y = yMax

        red = makeMaterial('Red',(1,0,0),(1,1,1),1)
        origin = (xMax,y,z)
        bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
        bpy.ops.transform.resize(value=(size, size, size))
        setMaterial(bpy.context.object, red)
        origin = (xMin,y,z)
        bpy.ops.mesh.primitive_uv_sphere_add(location=origin)
        bpy.ops.transform.resize(value=(size, size, size))
        setMaterial(bpy.context.object, red)
'''