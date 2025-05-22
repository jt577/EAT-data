import bpy
import json
import math

# Load atom data
with open('/Users/justint/Library/CloudStorage/OneDrive-Personal/Desktop/Academic Stuff/Arias Research/EAT/Paper_EAT_2025/Blender/EAT_noMag/run8+1_6/atom_data.json', 'r') as pkl_file:
    atom_data = json.load(pkl_file)

# Load element colors
with open('/Users/justint/Library/CloudStorage/OneDrive-Personal/Desktop/Academic Stuff/Arias Research/EAT/Paper_EAT_2025/Blender/EAT_noMag/run8+1_6/element_colors.json', 'r') as pkl_file:
    element_colors = json.load(pkl_file)

# Load element sizes (if needed)
with open('/Users/justint/Library/CloudStorage/OneDrive-Personal/Desktop/Academic Stuff/Arias Research/EAT/Paper_EAT_2025/Blender/EAT_noMag/run8+1_6/element_sizes.json', 'r') as pkl_file:
    element_sizes = json.load(pkl_file)

atom_objects = {}
sphere_scale = 0.75

# Access the first time step
initial_step = atom_data[0]

def create_atom_material(element_colors, atom_id):
    material = bpy.data.materials.new(name=f'AtomMaterial_{atom_id}')
    material.use_nodes = True
    nodes = material.node_tree.nodes
    links = material.node_tree.links
    nodes.clear()
    
    # Output node
    output_node = nodes.new(type='ShaderNodeOutputMaterial')
    output_node.location = (400, 0)
    
    # Principled BSDF shader
    shader_node = nodes.new(type='ShaderNodeBsdfPrincipled')
    shader_node.location = (200, 0)
    shader_node.inputs['Metallic'].default_value = 0.3
    shader_node.inputs['Roughness'].default_value = 0.25
    
    # Start with a zero vector (black color)
    base_color_node = nodes.new(type='ShaderNodeRGB')
    base_color_node.outputs[0].default_value = (0.0, 0.0, 0.0, 1.0)
    base_color_node.location = (-600, 0)
    
    previous_node = base_color_node
    
    # Create color and vector math nodes for each element
    for idx, (element, color) in enumerate(element_colors.items()):
        # RGB node for element color
        color_node = nodes.new(type='ShaderNodeRGB')
        color_node.outputs[0].default_value = (*color, 1.0)
        color_node.label = f'{element}_Color'
        color_node.name = f'{element}_Color'
        color_node.location = (-800, -200 * (idx + 1))
    
        # Value node for weight
        weight_node = nodes.new(type='ShaderNodeValue')
        weight_node.outputs[0].default_value = 0.0
        weight_node.label = f'{element}_Weight'
        weight_node.name = f'{element}_Weight'
        weight_node.location = (-1000, -200 * (idx + 1))
    
        # Vector Math node for scaling (multiply color by weight)
        scale_node = nodes.new(type='ShaderNodeVectorMath')
        scale_node.operation = 'SCALE'
        scale_node.location = (-600, -200 * (idx + 1))
    
        links.new(color_node.outputs[0], scale_node.inputs[0])
        links.new(weight_node.outputs[0], scale_node.inputs[3])
    
        # Vector Math node for adding scaled colors
        add_node = nodes.new(type='ShaderNodeVectorMath')
        add_node.operation = 'ADD'
        add_node.location = (-400, -100 * idx)
    
        links.new(previous_node.outputs[0], add_node.inputs[0])
        links.new(scale_node.outputs[0], add_node.inputs[1])
    
        previous_node = add_node
    
    # Connect to shader
    links.new(previous_node.outputs[0], shader_node.inputs['Base Color'])
    links.new(shader_node.outputs[0], output_node.inputs['Surface'])
    
    return material

# Create spheres and assign materials
for atom in initial_step['atoms']:
    atom_id = atom['id']
    position = atom['position']
    composition = atom['composition']

    # Add sphere
    bpy.ops.mesh.primitive_uv_sphere_add(radius=1.0, location=position)
    sphere = bpy.context.object
    sphere.name = f'Atom_{atom_id}'
    bpy.ops.object.shade_smooth()

    # Scale sphere based on dominant element
    if composition:
        dominant_element = max(composition, key=composition.get)
        radius = element_sizes.get(dominant_element, 1.0)
    else:
        radius = 1.0
    sphere.scale = (sphere_scale*radius,)*3

    # Create and assign material
    mat = create_atom_material(element_colors, atom_id)
    sphere.data.materials.append(mat)
    atom_objects[atom_id] = sphere

def update_atoms(time_step, frame_number):
    for atom in time_step['atoms']:
        atom_id = atom['id']
        position = atom['position']
        composition = atom['composition']
        sphere = atom_objects[atom_id]

        # Update position
        sphere.location = position
        sphere.keyframe_insert(data_path='location', frame=frame_number)

        # Update material weights
        mat = sphere.data.materials[0]
        nodes = mat.node_tree.nodes
        for element in element_colors.keys():
            weight = composition.get(element, 0.0)
            weight_node = nodes.get(f'{element}_Weight')
            if weight_node:
                weight_node.outputs[0].default_value = weight
                weight_node.outputs[0].keyframe_insert('default_value', frame=frame_number)

        # Update scale
        if composition:
            dominant_element = max(composition, key=composition.get)
            radius = element_sizes.get(dominant_element, 1.0)
        else:
            radius = 1.0
        sphere.scale = (sphere_scale*radius,)*3
        sphere.keyframe_insert(data_path='scale', frame=frame_number)

# Animate through all timesteps
for time_step in atom_data:
    frame = time_step['frame']
    bpy.context.scene.frame_set(frame)
    update_atoms(time_step, frame)

# Set scene frames, FPS, and resolution
scene = bpy.context.scene
scene.frame_start = atom_data[0]['frame']
scene.frame_end = atom_data[-1]['frame']
scene.render.fps = 6
scene.render.resolution_x = 1600
scene.render.resolution_y = 2600
scene.render.resolution_percentage = 100

# Add camera at specified position and rotation
bpy.ops.object.camera_add(
    location=(-11.1448, -28.7734, 10.021),
    rotation=(
        math.radians(72.7341),
        math.radians(0.001828),
        math.radians(-29.0343)
    )
)
camera = bpy.context.object
scene.camera = camera

# Add first point light
bpy.ops.object.light_add(type='POINT', location=(-7.31237, -21.731, 13.1237))
light1 = bpy.context.object
light1.data.energy = 30000
light1.data.shadow_soft_size = 10

# Add second point light
bpy.ops.object.light_add(type='POINT', location=(17.3717, 14.708, 0.328533))
light2 = bpy.context.object
light2.data.energy = 2000
light2.data.shadow_soft_size = 10
