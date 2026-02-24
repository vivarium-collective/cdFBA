from dataclasses import dataclass, field
from bigraph_schema.schema import Node, Map, List, Float
from bigraph_schema.methods import apply
from cdFBA.processes import register_processes
from plum import dispatch

#====================
#Volumetric Dataclass
#====================
@dataclass(kw_only=True)
class Volumetric(Node):
    concentrations: Map = field(default_factory= lambda: Map(_value=Float()))
    counts: Map = field(default_factory= lambda: Map(_value=Float()))
    volume: Float = field(default_factory= Float)

def conditional_apply(schema, current, update, key, path):
    """only applies update if specific key is found in the update"""
    if key in update:
        if isinstance(schema, Node):
            subschema = getattr(schema, key)
        else:
            subschema = schema[key]
        applied, merges = apply(
            subschema,
            current[key],
            update[key],
            path + (key, ),
        )
    else:
        applied = current[key]
        merges = []

    return applied, merges



def volumetric_update(schema, current, update, path):
    """applies update to a volumetric type"""
    if not "concentration" in update:
        updated_counts, merges_counts = conditional_apply(schema, current, update, 'counts', path)
        updated_volume, merges_volume = conditional_apply(schema, current, update, 'volume', path)
        updated_concentrations = {}
        merges = merges_volume + merges_counts
        for key, counts in updated_counts.items():
            updated_concentrations[key] = counts / updated_volume

    else:
        if "volume" in update:
            raise ValueError("Cannot apply volume and concentration updates at the same time")
        if "counts" in update:
            raise ValueError("Cannot apply count and concentration updates at the same time")

        updated_concentrations, merges = conditional_apply(schema, current, update, 'concentrations', path)
        updated_volume = current["volume"]
        updated_counts = {}
        for key, concentration in updated_concentrations.items():
            updated_counts[key] =  updated_volume * concentration

    applied = {
        'counts': updated_counts,
        'concentrations': updated_concentrations,
        'volume': updated_volume,
    }

    return applied, merges

@apply.dispatch
def apply(schema: Volumetric, current, update, path):
    return volumetric_update(schema, current, update, path)

# @dispatch
# def resolve(current: Volumetric, update: Map, path=None):


bounds_type = {
    "lower": "maybe[float]",
    "upper": "maybe[float]"
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
    "bounds": "map[bounds]",
    "kinetics": "map",
}

threshold_type = {
    "type": "string",  # add or remove
    "substrate": "string",  # substrate or species to monitor
    "range": {
        "upper": "maybe[float]",
        "lower": "maybe[float]",
    },
    "parent": "string",  # bane of parent species
    "name": "string",
    "changes": "dfba_changes",
    "mass": "float"
}


def register_types(core):
    core.register_type("bounds", bounds_type)
    core.register_type("volumetric", Volumetric)
    core.register_type("threshold", threshold_type)
    core.register_type("dfba_changes", dfba_changes_type)

    return register_processes(core)
