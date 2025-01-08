from process_bigraph import ProcessTypes
from cdFBA.processes import register_processes

def apply_non_negative(schema, current, update, core):
    new_value = current + update
    return max(0, new_value)


def set_update(schema, current, update, core):
    return update


positive_float = {
    '_type': 'positive_float',
    '_inherit': 'float',
    '_apply': apply_non_negative}

set_float = {
    '_type': 'set_float',
    '_inherit': 'float',
    '_apply': set_update
}

bounds_type = {
    'lower': 'maybe[float]',
    'upper': 'maybe[float]'}


particle_type = {
    'id': 'string',
    'position': 'tuple[float,float]',
    'size': 'float',
    'local': 'map[float]',
    'exchange': 'map[float]',    # {mol_id: delta_value}
}

def register_types(core):
    core.register('positive_float', positive_float)
    core.register('set_float', set_float)
    core.register('bounds', bounds_type)
    core.register('particle', particle_type)
    
    return register_processes(core)