from cdFBA.processes import register_processes

def apply_non_negative(schema, current, update, core):
    new_value = current + update
    return max(0, new_value)

def set_update(schema, current, update, top_schema, top_state, path, core):
    return update

def conditional_apply(schema, current, update, key, core):
    if key in update:
        applied = core.apply(
            schema[key],
            current[key],
            update[key])
    else:
        applied = current[key]

    return applied

def volumetric_update(schema, current, update, top_schema, top_state, path, core):
    updated_counts = conditional_apply(schema, current, update, 'counts', core)
    updated_volume = conditional_apply(schema, current, update, 'volume', core)
    updated_concentrations = {}

    for key, counts in updated_counts.items():
        updated_concentrations[key] = counts / updated_volume

    applied = {
        'counts': updated_counts,
        'concentrations': updated_concentrations,
        'volume': updated_volume,
    }
    
    return applied

positive_float = {
    "_type": "positive_float",
    "_inherit": "float",
    "_apply": apply_non_negative}

set_float = {
    "_type": "set_float",
    "_inherit": "float",
    "_apply": set_update
}

bounds_type = {
    "lower": "maybe[float]",
    "upper": "maybe[float]"}


particle_type = {
    "id": "string",
    "position": "tuple[float,float]",
    "size": "float",
    "local": "map[float]",
    "exchange": "map[float]",    # {mol_id: delta_value}
}

volumetric_type = {
    "concentrations":"map[float]",
    "counts":"map[float]",
    "volume": {"_type": "float", "_default": 1.0},
    "_apply": volumetric_update
}

chemostat_type = {
    "concentration": "float"
}

dfba_launch_type = {
    "model_file": "string",
    "config": "map",
    "biomass": "float",
    "parent": "string"
}

dfba_changes_type = {
    "gene_knockout": "maybe[list]",
    "reaction_knockout": "maybe[list]",
    "bounds": "any",
    "kinetics": "any",
}

threshold_type = {
    "type": "string", #add or remove
    "substrate": "string", #substrate or species to monitor
    "range": {
        "upper": "maybe[float]",
        "lower": "maybe[float]",
    },
    "parent": "string", #bane of parent species
    "name": "string",
    "changes": "dfba_changes",
    "mass": "float"
}

def register_types(core):
    core.register("positive_float", positive_float)
    core.register("set_float", set_float)
    core.register("bounds", bounds_type)
    core.register("particle", particle_type)
    core.register("volumetric", volumetric_type)
    core.register("threshold", threshold_type)
    core.register("dfba_changes", dfba_changes_type)
    
    return register_processes(core)
