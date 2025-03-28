from cdFBA.processes import register_processes
#TODO -- add global store name variables

def apply_non_negative(schema, current, update, core):
    new_value = current + update
    return max(0, new_value)

def set_update(schema, current, update, core):
    return update

def volumetric_update(schema, current, update, core):
    counts = update.get("counts")
    volume = current.get("volume")

    new_counts = current["counts"]
    new_concentrations = current["concentrations"]
    if counts:
        for key, value in counts.items():
            new = new_counts[key] + value
            new_counts[key] = new
            new_concentrations[key] = new/volume

    return {
        "counts": new_counts,
        "concentrations": new_concentrations,
        "volume": volume,
    }

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
    "gene_knockout": "maybe[list[string]]",
    "reaction_knockout": "maybe[list[string]]",
}

threshold_type = {
    "substrate": "string",
    "range": {
        "upper": "maybe[float]",
        "lower": "maybe[float]",
    },
    "parent": "string",
    "name": "string",
    "changes": "maybe[dfba_changes]"
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