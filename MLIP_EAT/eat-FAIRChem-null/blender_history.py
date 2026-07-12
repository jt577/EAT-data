import argparse
import json
import sys
from pathlib import Path

import bpy
from mathutils import Vector


DEFAULT_HISTORY_PATH = Path("/Users/justint/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Arias Research/EAT/2026/04.2026/CuPd_nulls/all_history_target=Cu0.750_Pd0.250_i=0.json")
DEFAULT_FRAMES_PER_STATE = 2
DEFAULT_FPS = 40
DEFAULT_RESOLUTION_X = 1000
DEFAULT_RESOLUTION_Y = 1000
DEFAULT_ATOM_RADIUS = 0.55
DEFAULT_CELL_LINE_WIDTH = 0.05
DEFAULT_LIGHTING_STRENGTH = 0.1
DEFAULT_COLOR_WEIGHT_POWER = 1.0
COLLECTION_NAME = "EATHistory"
CELL_EDGE_INDICES = (
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 4),
    (1, 5),
    (2, 4),
    (2, 6),
    (3, 5),
    (3, 6),
    (4, 7),
    (5, 7),
    (6, 7),
)
DEFAULT_SPECIES_COLORS = {
    "Na": (0.96, 0.74, 0.10, 1.0),
    "Cl": (0.12, 0.76, 0.20, 1.0),
    "Li": (1.0, 0.25, 0.20, 1.0),
    "O": (1.0, 0.22, 0.28, 1.0),
    "F": (0.12, 0.98, 0.82, 1.0),
    "S": (1.0, 0.92, 0.10, 1.0),
    "Mg": (0.18, 0.78, 1.0, 1.0),
    "Al": (0.88, 0.88, 0.92, 1.0),
    "Si": (0.72, 0.48, 1.0, 1.0),
    "P": (1.0, 0.55, 0.08, 1.0),
    "K": (0.74, 0.32, 1.0, 1.0),
    "Ca": (0.32, 0.96, 0.84, 1.0),
}
FALLBACK_PALETTE = (
    (1.0, 0.22, 0.28, 1.0),
    (0.14, 0.86, 1.0, 1.0),
    (0.96, 0.74, 0.10, 1.0),
    (0.12, 0.76, 0.20, 1.0),
    (0.82, 0.35, 1.0, 1.0),
    (1.0, 0.56, 0.10, 1.0),
)

ANIMATION_CACHE = {}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--history",
        default=str(DEFAULT_HISTORY_PATH),
        help="Path to all_history_i=*.json written by the EAT optimizer.",
    )
    parser.add_argument(
        "--frames-per-state",
        type=int,
        default=DEFAULT_FRAMES_PER_STATE,
        help="How many rendered frames each optimizer state occupies.",
    )
    parser.add_argument("--fps", type=int, default=DEFAULT_FPS)
    parser.add_argument("--resolution-x", type=int, default=DEFAULT_RESOLUTION_X)
    parser.add_argument("--resolution-y", type=int, default=DEFAULT_RESOLUTION_Y)
    parser.add_argument(
        "--atom-radius",
        type=float,
        default=DEFAULT_ATOM_RADIUS,
        help="Base radius in angstrom-like Blender units for each atom sphere.",
    )
    parser.add_argument(
        "--cell-line-width",
        type=float,
        default=DEFAULT_CELL_LINE_WIDTH,
        help="Bevel depth for the unit-cell edges.",
    )
    parser.add_argument(
        "--lighting-strength",
        type=float,
        default=DEFAULT_LIGHTING_STRENGTH,
        help="Global multiplier for the full lighting rig around the cell.",
    )
    parser.add_argument(
        "--color-weight-power",
        type=float,
        default=DEFAULT_COLOR_WEIGHT_POWER,
        help="Nonlinear power used to emphasize dominant species colors before renormalizing.",
    )
    if "--" in sys.argv:
        cli_args = sys.argv[sys.argv.index("--") + 1 :]
    else:
        cli_args = []
    return parser.parse_args(cli_args)


def load_history(path):
    with Path(path).expanduser().open("r", encoding="utf-8") as handle:
        return json.load(handle)


def set_node_input(node, candidate_names, value):
    for name in candidate_names:
        socket = node.inputs.get(name)
        if socket is not None:
            socket.default_value = value
            return


def build_species_palette(species_order):
    palette = {}
    for idx, species in enumerate(species_order):
        if species in DEFAULT_SPECIES_COLORS:
            palette[species] = DEFAULT_SPECIES_COLORS[species]
        else:
            palette[species] = FALLBACK_PALETTE[idx % len(FALLBACK_PALETTE)]
    return palette


def snapshot_signature(snapshot):
    return json.dumps(
        {
            "iter": snapshot.get("iter"),
            "cell": snapshot.get("lattice_vectors_angstrom"),
            "positions": snapshot.get("cart_positions_angstrom"),
            "weights": snapshot.get("species_weights"),
            "occupancies": snapshot.get("site_occupancies"),
        },
        sort_keys=True,
    )


def build_snapshot_from_parallel(history, index):
    return {
        "iter": history.get("iter", [index + 1])[index],
        "frac_coords": history["frac_coords"][index],
        "cart_positions_angstrom": history["cart_positions_angstrom"][index],
        "lattice_vectors_angstrom": history["lattice_vectors_angstrom"][index],
        "species_weights": history["species_weights"][index],
        "site_occupancies": history.get("site_occupancies", [None] * len(history["species_weights"]))[index],
        "dominant_species": history.get("dominant_species", [None] * len(history["species_weights"]))[index],
        "cost": history.get("cost", [None] * len(history["species_weights"]))[index],
        "total_energy": history.get("total_energy", [None] * len(history["species_weights"]))[index],
        "formation_energy": history.get("formation_energy", [None] * len(history["species_weights"]))[index],
        "energy_above_hull": history.get("energy_above_hull", [None] * len(history["species_weights"]))[index],
        "entropy_per_atom": history.get("entropy_per_atom", [None] * len(history["species_weights"]))[index],
        "cell_volume": history.get("cell_volume", [None] * len(history["species_weights"]))[index],
    }


def append_state(frames, snapshot):
    if snapshot is None:
        return
    signature = snapshot_signature(snapshot)
    if frames and frames[-1]["signature"] == signature:
        return
    frames.append({"snapshot": snapshot, "signature": signature})


def append_run(frames, run_record, include_initial=True):
    if not isinstance(run_record, dict):
        return
    if include_initial:
        append_state(frames, run_record.get("initial_state"))
    for snapshot in run_record.get("trajectory", []):
        append_state(frames, snapshot)


def find_accepted_trial(escape_round):
    best_trial = escape_round.get("best_trial")
    for trial in escape_round.get("trial_records", []):
        if trial.get("accepted"):
            return trial
    if best_trial is not None:
        for trial in escape_round.get("trial_records", []):
            if trial.get("trial_index") == best_trial:
                return trial
    return None


def build_frames(history):
    frames = []
    raw_stages = history.get("stages", [])
    if raw_stages:
        append_state(frames, history.get("initial_state"))
        for stage in raw_stages:
            runs = list(stage.get("optimization_runs", []))
            escape_rounds = list(stage.get("escape_rounds", []))

            # In the current null syntropizer, stage runs are recorded in chronological order
            # before any accepted escape is applied. So we should play each recorded run once,
            # then append the accepted shock/relaxation sequence, rather than replaying an older
            # pre-shock run after the escape.
            for run_index, run in enumerate(runs):
                append_run(frames, run, include_initial=(run_index == 0))

            for escape_round in escape_rounds:
                accepted_trial = find_accepted_trial(escape_round) if escape_round.get("accepted") else None
                if accepted_trial is not None:
                    append_state(frames, accepted_trial.get("shock_state"))
                    append_run(frames, accepted_trial.get("relaxation", {}), include_initial=False)

        append_run(frames, history.get("final_relaxation", {}), include_initial=True)
    else:
        raw_trajectory = history.get("trajectory", [])
        if raw_trajectory:
            snapshots = raw_trajectory
        else:
            n_steps = len(history.get("species_weights", []))
            snapshots = [build_snapshot_from_parallel(history, idx) for idx in range(n_steps)]

        previous_signature = None
        for snapshot in snapshots:
            signature = snapshot_signature(snapshot)
            if signature == previous_signature:
                continue
            frames.append({"snapshot": snapshot, "signature": signature})
            previous_signature = signature

    for index, frame_record in enumerate(frames):
        frame_record["cumulative_iteration"] = index
    return frames


def clear_scene():
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for mesh in list(bpy.data.meshes):
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)
    for material in list(bpy.data.materials):
        if material.users == 0:
            bpy.data.materials.remove(material)
    for curve in list(bpy.data.curves):
        if curve.users == 0:
            bpy.data.curves.remove(curve)
    for camera in list(bpy.data.cameras):
        if camera.users == 0:
            bpy.data.cameras.remove(camera)
    for light in list(bpy.data.lights):
        if light.users == 0:
            bpy.data.lights.remove(light)
    for collection in list(bpy.data.collections):
        if collection.users == 0:
            bpy.data.collections.remove(collection)


def ensure_collection(name):
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    return collection


def move_object_to_collection(obj, target_collection):
    if target_collection not in obj.users_collection:
        target_collection.objects.link(obj)
    for existing_collection in list(obj.users_collection):
        if existing_collection != target_collection:
            existing_collection.objects.unlink(obj)


def create_atom_material(name, rgba):
    material = bpy.data.materials.new(name=name)
    material.use_nodes = True
    if hasattr(material, "blend_method"):
        material.blend_method = "BLEND"
    if hasattr(material, "shadow_method"):
        material.shadow_method = "HASHED"

    nodes = material.node_tree.nodes
    links = material.node_tree.links
    nodes.clear()

    output_node = nodes.new(type="ShaderNodeOutputMaterial")
    output_node.location = (260, 0)

    shader_node = nodes.new(type="ShaderNodeBsdfPrincipled")
    shader_node.location = (20, 0)
    set_node_input(shader_node, ("Base Color",), rgba)
    set_node_input(shader_node, ("Alpha",), rgba[3])
    set_node_input(shader_node, ("Metallic",), 0.0)
    set_node_input(shader_node, ("Roughness",), 0.20)
    set_node_input(shader_node, ("Specular IOR Level", "Specular"), 0.42)
    set_node_input(shader_node, ("Coat Weight", "Clearcoat"), 0.12)
    set_node_input(shader_node, ("Coat Roughness", "Clearcoat Roughness"), 0.12)

    links.new(shader_node.outputs[0], output_node.inputs["Surface"])
    return material


def create_color_material(name, rgba):
    material = bpy.data.materials.new(name=name)
    material.use_nodes = True
    nodes = material.node_tree.nodes
    links = material.node_tree.links
    nodes.clear()

    output_node = nodes.new(type="ShaderNodeOutputMaterial")
    shader_node = nodes.new(type="ShaderNodeBsdfPrincipled")
    set_node_input(shader_node, ("Base Color",), rgba)
    set_node_input(shader_node, ("Emission Color", "Emission"), (*rgba[:3], 1.0))
    set_node_input(shader_node, ("Emission Strength",), 0.0)
    set_node_input(shader_node, ("Roughness",), 0.45)

    links.new(shader_node.outputs[0], output_node.inputs["Surface"])
    return material


def create_atom_objects(collection, atom_count, initial_snapshot, atom_radius, palette, species_order):
    atom_objects = []
    dominant_species = initial_snapshot.get("dominant_species", [])
    weights = initial_snapshot["species_weights"]
    occupancies = initial_snapshot.get("site_occupancies")
    for atom_index in range(atom_count):
        rgba = blended_rgba(weights[atom_index], species_order, palette)
        occupancy = 1.0
        if occupancies is not None and atom_index < len(occupancies):
            occupancy = max(0.0, min(1.0, float(occupancies[atom_index])))
        rgba = (rgba[0], rgba[1], rgba[2], occupancy)
        bpy.ops.mesh.primitive_uv_sphere_add(
            radius=1.0,
            location=(0.0, 0.0, 0.0),
            segments=48,
            ring_count=24,
        )
        sphere = bpy.context.object
        bpy.ops.object.shade_smooth()
        sphere.scale = (atom_radius, atom_radius, atom_radius)
        suffix = dominant_species[atom_index] if atom_index < len(dominant_species) else "mix"
        sphere.name = f"Atom_{atom_index:03d}_{suffix}"
        material = create_atom_material(f"AtomMaterial_{atom_index:03d}", rgba)
        sphere.data.materials.clear()
        sphere.data.materials.append(material)
        move_object_to_collection(sphere, collection)
        atom_objects.append(sphere)
    return atom_objects


def create_cell_edge_objects(collection, line_width):
    edge_material = create_color_material("CellEdgeMaterial", (1.0, 1.0, 1.0, 1.0))
    edge_objects = []
    for edge_index, _edge in enumerate(CELL_EDGE_INDICES):
        curve = bpy.data.curves.new(f"CellEdgeCurve_{edge_index:02d}", type="CURVE")
        curve.dimensions = "3D"
        curve.bevel_depth = line_width
        curve.bevel_resolution = 3
        spline = curve.splines.new("POLY")
        spline.points.add(1)
        spline.points[0].co = (0.0, 0.0, 0.0, 1.0)
        spline.points[1].co = (0.0, 0.0, 0.0, 1.0)
        edge_object = bpy.data.objects.new(f"CellEdge_{edge_index:02d}", curve)
        edge_object.data.materials.append(edge_material)
        collection.objects.link(edge_object)
        edge_objects.append(edge_object)
    return edge_objects


def compute_cell_corners(cell_vectors):
    a_vec = Vector(cell_vectors[0])
    b_vec = Vector(cell_vectors[1])
    c_vec = Vector(cell_vectors[2])
    origin = Vector((0.0, 0.0, 0.0))
    return (
        origin,
        a_vec,
        b_vec,
        c_vec,
        a_vec + b_vec,
        a_vec + c_vec,
        b_vec + c_vec,
        a_vec + b_vec + c_vec,
    )


def update_cell_edges(edge_objects, cell_vectors, frame_number=None, keyframe=False):
    corners = compute_cell_corners(cell_vectors)
    for edge_object, (start_index, end_index) in zip(edge_objects, CELL_EDGE_INDICES):
        spline = edge_object.data.splines[0]
        spline.points[0].co = (*corners[start_index], 1.0)
        spline.points[1].co = (*corners[end_index], 1.0)
        if keyframe and frame_number is not None:
            spline.points[0].keyframe_insert(data_path="co", frame=frame_number)
            spline.points[1].keyframe_insert(data_path="co", frame=frame_number)
        if hasattr(edge_object.data, "update_tag"):
            edge_object.data.update_tag()


def blended_rgba(weights, species_order, palette, power=DEFAULT_COLOR_WEIGHT_POWER):
    total = sum(max(0.0, float(weight)) for weight in weights)
    if total <= 1.0e-12:
        return (0.55, 0.55, 0.58, 1.0)

    normalized_weights = [max(0.0, float(weight)) / total for weight in weights]
    nonlinear_weights = []
    power = max(1.0e-6, float(power))
    for weight in normalized_weights:
        nonlinear_weights.append(weight ** power)

    nonlinear_total = sum(nonlinear_weights)
    if nonlinear_total <= 1.0e-12:
        nonlinear_weights = [max(0.0, float(weight)) for weight in weights]
        nonlinear_total = sum(nonlinear_weights)

    rgb = [0.0, 0.0, 0.0]
    for weight, species in zip(nonlinear_weights, species_order):
        color = palette[species]
        rgb[0] += weight * color[0]
        rgb[1] += weight * color[1]
        rgb[2] += weight * color[2]

    rgb = [channel / nonlinear_total for channel in rgb]
    return (rgb[0], rgb[1], rgb[2], 1.0)


def apply_frame_record(frame_record, frame_number=None, keyframe=False):
    snapshot = frame_record["snapshot"]
    positions = snapshot["cart_positions_angstrom"]
    weights = snapshot["species_weights"]
    occupancies = snapshot.get("site_occupancies")
    cell_vectors = snapshot["lattice_vectors_angstrom"]
    species_order = ANIMATION_CACHE["species_order"]
    palette = ANIMATION_CACHE["palette"]
    atom_radius = ANIMATION_CACHE["atom_radius"]
    color_weight_power = ANIMATION_CACHE["color_weight_power"]

    for atom_index, (atom_object, position, atom_weights) in enumerate(
        zip(ANIMATION_CACHE["atom_objects"], positions, weights)
    ):
        atom_object.location = position
        rgba = blended_rgba(atom_weights, species_order, palette, color_weight_power)
        total_weight = sum(max(0.0, float(weight)) for weight in atom_weights)
        occupancy = None
        if occupancies is not None and atom_index < len(occupancies):
            occupancy = max(0.0, min(1.0, float(occupancies[atom_index])))
        if occupancy is None:
            occupancy = max(0.0, min(1.0, float(total_weight)))
        if occupancy <= 1.0e-12:
            atom_object.scale = (0.0, 0.0, 0.0)
            atom_object.hide_render = True
            atom_object.hide_viewport = True
            if keyframe and frame_number is not None:
                atom_object.keyframe_insert(data_path="location", frame=frame_number)
                atom_object.keyframe_insert(data_path="scale", frame=frame_number)
                atom_object.keyframe_insert(data_path="hide_render", frame=frame_number)
                atom_object.keyframe_insert(data_path="hide_viewport", frame=frame_number)
            continue
        atom_object.scale = (atom_radius, atom_radius, atom_radius)
        material = atom_object.active_material
        principled = next(
            node for node in material.node_tree.nodes if node.type == "BSDF_PRINCIPLED"
        )
        display_rgba = (rgba[0], rgba[1], rgba[2], occupancy)
        set_node_input(principled, ("Base Color",), display_rgba)
        set_node_input(principled, ("Alpha",), occupancy)
        atom_object.hide_render = False
        atom_object.hide_viewport = False
        if keyframe and frame_number is not None:
            atom_object.keyframe_insert(data_path="location", frame=frame_number)
            atom_object.keyframe_insert(data_path="scale", frame=frame_number)
            atom_object.keyframe_insert(data_path="hide_render", frame=frame_number)
            atom_object.keyframe_insert(data_path="hide_viewport", frame=frame_number)
            principled.inputs["Base Color"].keyframe_insert(
                data_path="default_value",
                frame=frame_number,
            )
            alpha_socket = principled.inputs.get("Alpha")
            if alpha_socket is not None:
                alpha_socket.keyframe_insert(data_path="default_value", frame=frame_number)

    update_cell_edges(
        ANIMATION_CACHE["cell_edges"],
        cell_vectors,
        frame_number=frame_number,
        keyframe=keyframe,
    )

def remove_existing_handlers():
    for handler in list(bpy.app.handlers.frame_change_pre):
        if getattr(handler, "_eat_history_handler", False):
            bpy.app.handlers.frame_change_pre.remove(handler)

def bake_animation():
    frames_per_state = ANIMATION_CACHE["frames_per_state"]
    for state_index, frame_record in enumerate(ANIMATION_CACHE["frames"]):
        frame_start = 1 + state_index * frames_per_state
        apply_frame_record(frame_record, frame_number=frame_start, keyframe=True)

        if frames_per_state > 1:
            hold_frame = frame_start + frames_per_state - 1
            apply_frame_record(frame_record, frame_number=hold_frame, keyframe=True)


def compute_view_setup(initial_snapshot):
    cell_vectors = initial_snapshot["lattice_vectors_angstrom"]
    corners = compute_cell_corners(cell_vectors)
    center = sum(corners, Vector((0.0, 0.0, 0.0))) / len(corners)
    extents = [max(corner[axis] for corner in corners) for axis in range(3)]
    mins = [min(corner[axis] for corner in corners) for axis in range(3)]
    span = max(extents[axis] - mins[axis] for axis in range(3))
    span = max(span, 1.0)
    camera_location = center + Vector((-1.85 * span, -2.15 * span, 1.18 * span))
    return center, camera_location, span


def setup_world():
    world = bpy.context.scene.world
    world.use_nodes = True
    background = world.node_tree.nodes["Background"]
    background.inputs[0].default_value = (0.0, 0.0, 0.0, 1.0)
    background.inputs[1].default_value = 0.0


def create_camera_and_lights(collection, target_location, camera_location, span, lighting_strength):
    bpy.ops.object.camera_add(location=tuple(camera_location))
    camera = bpy.context.object
    camera.name = "HistoryCamera"
    move_object_to_collection(camera, collection)

    target = bpy.data.objects.new("CameraTarget", None)
    target.location = target_location
    collection.objects.link(target)

    constraint = camera.constraints.new(type="TRACK_TO")
    constraint.target = target
    constraint.track_axis = "TRACK_NEGATIVE_Z"
    constraint.up_axis = "UP_Y"

    strength = max(0.0, float(lighting_strength))
    # Symmetric lighting shell: same angular coverage from every direction,
    # with equal-energy lights on upper, equatorial, and lower rings.
    ring_radius = 1.65 * span
    z_levels = (
        ("Upper", 0.78, ((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0))),
        ("Mid", 0.0, ((0.707, 0.707), (-0.707, 0.707), (-0.707, -0.707), (0.707, -0.707))),
        ("Lower", -0.78, ((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0))),
    )
    light_specs = []
    for ring_name, z_dir, xy_dirs in z_levels:
        for idx, (x_dir, y_dir) in enumerate(xy_dirs):
            direction = Vector((x_dir, y_dir, z_dir)).normalized()
            light_specs.append(
                (
                    f"{ring_name}Light_{idx:02d}",
                    direction * ring_radius,
                    2800.0,
                    1.20 * span,
                    1.20 * span,
                )
            )

    for name, offset, energy, size, size_y in light_specs:
        bpy.ops.object.light_add(type="AREA", location=tuple(target_location + offset))
        light = bpy.context.object
        light.name = name
        light.data.energy = energy * strength
        light.data.shape = "DISK" if abs(size - size_y) < 1.0e-8 else "ELLIPSE"
        light.data.size = size
        if hasattr(light.data, "size_y"):
            light.data.size_y = size_y
        if hasattr(light.data, "spread"):
            light.data.spread = 3.14159

        track = light.constraints.new(type="TRACK_TO")
        track.target = target
        track.track_axis = "TRACK_NEGATIVE_Z"
        track.up_axis = "UP_Y"
        move_object_to_collection(light, collection)

    return camera


def configure_scene(args, frame_count):
    scene = bpy.context.scene
    scene.frame_start = 1
    scene.frame_end = max(1, args.frames_per_state * (frame_count - 1) + 1)
    scene.render.fps = args.fps
    scene.render.resolution_x = args.resolution_x
    scene.render.resolution_y = args.resolution_y
    scene.render.resolution_percentage = 100
    for engine_name in ("BLENDER_EEVEE_NEXT", "BLENDER_EEVEE", "CYCLES"):
        try:
            scene.render.engine = engine_name
            break
        except TypeError:
            continue
    if hasattr(scene, "eevee"):
        scene.eevee.taa_render_samples = 128
        scene.eevee.taa_samples = 64
    if hasattr(scene, "cycles"):
        scene.cycles.samples = 128
        scene.cycles.preview_samples = 32
    if hasattr(scene.render, "film_transparent"):
        scene.render.film_transparent = False
    setup_world()


def main():
    args = parse_args()
    history = load_history(args.history)
    species_order = history.get("species_order") or history.get("selected_elements")
    if not species_order:
        raise RuntimeError("History JSON is missing 'species_order' and 'selected_elements'.")

    frames = build_frames(history)
    if not frames:
        raise RuntimeError("No animation frames could be constructed from the history JSON.")

    clear_scene()
    collection = ensure_collection(COLLECTION_NAME)
    initial_snapshot = frames[0]["snapshot"]
    center, camera_location, span = compute_view_setup(initial_snapshot)
    palette = build_species_palette(species_order)

    camera = create_camera_and_lights(
        collection,
        center,
        camera_location,
        span,
        args.lighting_strength,
    )
    bpy.context.scene.camera = camera

    atom_count = len(initial_snapshot["cart_positions_angstrom"])
    atom_objects = create_atom_objects(
        collection=collection,
        atom_count=atom_count,
        initial_snapshot=initial_snapshot,
        atom_radius=args.atom_radius,
        palette=palette,
        species_order=species_order,
    )
    cell_edges = create_cell_edge_objects(collection, args.cell_line_width)

    ANIMATION_CACHE.clear()
    ANIMATION_CACHE.update(
        {
            "frames": frames,
            "frames_per_state": max(1, args.frames_per_state),
            "atom_objects": atom_objects,
            "cell_edges": cell_edges,
            "species_order": species_order,
            "palette": palette,
            "atom_radius": args.atom_radius,
            "color_weight_power": args.color_weight_power,
        }
    )

    remove_existing_handlers()
    configure_scene(args, len(frames))
    bake_animation()
    bpy.context.scene.frame_set(1)
    apply_frame_record(frames[0])


if __name__ == "__main__":
    main()
