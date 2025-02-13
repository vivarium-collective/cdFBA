"""This module contains some methods to obtain the reaction map, initial conditions, and kinetic parameters needed
for dFBA simulations from the minimal medium requirements of the wild-type species.

CAUTION: The initial conditions, and kinetics dataframes provide default parameter values and must be changed as needed
CAUTION: Substrate names are different in BiGG and AGORA databases. These functions will not work with two models form
         different sources
TODO: add user input for substrates
"""
from copy import deepcopy

#TODO: make notebook demo for functions and make them pure functions

from cobra.io import load_model, read_sbml_model, load_json_model, load_yaml_model, load_matlab_model
from cobra.medium import minimal_medium
import pprint
import re
import copy

def make_cdfba_composite(model_dict, medium='exchange'):
    """Construct a cdfba composite spec with all exhange metabolites included.
    Parameters:
    model_dict : dict, dictionary with cdfba process names as keys and model name/path as values
    medium : str/list, if str, pick one of:
        'default' uses the default cobra model medium
        'minimal' uses the minimal medium for the model
        'exchange' uses all exchange fluxes for the model
    if list:
        a list of exchange reaction ids

    Returns:
    spec : cdfba composite spec
    """
    if isinstance(medium, str):

        models = {model_name: model_from_file(model) for model_name, model in model_dict.items()}

    pass

def model_from_file(model_file='textbook'):
    """Returns a cobra model from a model file path or BiGG Model ID
    Parameters:
        model_file: string, file path or BiGG Model ID
    Returns:
        model: cobra model
    """
    if ".xml" in model_file:
        model = read_sbml_model(model_file)
    elif ".json" in model_file:
        model = load_json_model(model_file)
    elif ".yaml" in model_file:
        model = load_yaml_model(model_file)
    elif ".mat" in model_file:
        model = load_matlab_model(model_file)
    elif isinstance(model_file, str):
        model = load_model(model_file)
    else:
        # error handling
        raise ValueError("Invalid model file")
    return model

def get_exchanges(model_file='textbook', medium_type='exchange'):
    model=model_from_file(model_file)
    if medium_type == 'default':
        medium = model.medium
    if medium_type == 'minimal':
        medium = minimal_medium(model, model.slim_optimize()).to_dict()
    if medium_type == 'exchange':
        medium = {reaction.id: reaction.upper_bound for reaction in model.exchanges}
        medium.update(model.medium)
        medium = medium
    return list(medium.keys())

def get_substrates(model_file='textbook', exchanges=None):
    """Returns a list of substrates from the model.
    Parameters:
    model : cobrapy model
    exchanges : list, a list of exchange reaction ids
    Returns:
    substrates : list, list of names of substrates required by the model organism
    """
    if exchanges is None:
        exchanges = get_exchanges(model_file)
    model = model_from_file(model_file)
    substrates = []
    for item in [getattr(model.reactions, i).name for i in exchanges]:
        match = re.match(r"(.*) exchange|exchange reaction for (.*)|Exchange of (.*)|echange reaction for (.*)", item, re.IGNORECASE)
        if match:
            substrates.append(match.group(1) or match.group(2) or match.group(3) or match.group(4))
        else:
            substrates.append(item)
    return substrates
        
def get_reaction_map(model_file='textbook', exchanges=None):
    """Returns a reaction_name_map dictionary from a medium dictionary as obtained
    from model.medium or cobra.medium.minimum_medium()
    Parameters:
        model :
        exchanges : list, list of names of substrates required by the model organism
    Returns:
        reaction_name_map : dict, maps substrate names to reactions
    """
    if exchanges is None:
        exchanges = get_exchanges(model_file)
    model = model_from_file(model_file)
    substrates = get_substrates(model_file, exchanges)
    ids = exchanges
    reaction_name_map = {}
    for i in range(len(substrates)):
        reaction_name_map[substrates[i]] = ids[i]
    return reaction_name_map
    
def get_kinetics(model_file='textbook', exchanges=None):
    """Returns default kinetic parameters dictionary. Values are tuples of the form (km, vmax)"""
    if exchanges is None:
        exchanges = get_exchanges(model_file)
    model = model_from_file(model_file)
    kinetics = {key: (0.5, 2.0) for key in get_substrates(model_file, exchanges)}
    return kinetics

def get_bounds(reaction_map):
    """Return dict of upper and lower bounds for each substrate exchange reaction"""
    return {reaction_map[key]: {'lower': -1000, 'upper': 1000} for key in reaction_map.keys()}

def get_objective_reaction(model_file = 'textbook'):
    """get a string with the name of the objective function of a cobra model
    Parameters:
        model: cobrapy model
    Returns:
        objective_reaction: string, name of the objective reaction (biomass reaction by default)
    """
    model = model_from_file(model_file)
    expression = f"{model.objective.expression}"
    match = re.search(r'1\.0\*([^\s]+)', expression)

    if match:
        objective_reaction = match.group(1)

    return objective_reaction

def get_initial_counts(model_dict, biomass=0.1, initial_value=20, exchanges=None):
    """Returns an initial condition dict based on medium
    Parameters:
        model: string, cobrapy model name
        substrates : list, list of names of substrates required by the model organism
        biomass : float, initial biomass for all species
        factor : float, factor to multiply minimum medium concentrations
    Returns:
        conditions : dict, initial conditions dictionary
    """
    all_substrates = []
    for model_name in model_dict.keys():
        model_file=model_dict[model_name]
        substrates = get_substrates(model_file, exchanges)
    conditions = {substrate:initial_value for substrate in substrates}
    biomasses = {model:biomass for model in model_dict.keys()}
    conditions = conditions | biomasses
    return conditions

def initial_environment(volume=1, initial_counts=None, species_list=None):
    """Construct initial shared environment store
    Parameters:
        volume : float, volume of the environment
        initial_counts : dict, initial counts of each substrate and species biomass in the environment
        species_list : list of strings, list of dfba species names (DFBA.config["name"])
    Returns:
        initial shared environment store spec
    """
    if initial_counts is None:
        if species_list is None:
            raise ValueError("Error: Please provide initial_counts or species_list")
        initial_counts = {
            "glucose": 80,
            "acetate": 0,
        }
        for species in species_list:
            initial_counts[species] = 0.5

    initial_concentration = {key:(count/volume) for key, count in initial_counts.items()}

    return {
        "volume": volume,
        "counts": initial_counts,
        "concentrations": initial_concentration
    }

def dfba_config(
        model_file="textbook",
        name=None,
        kinetics=None,
        reaction_map=None,
        biomass_identifier=None,
        bounds=None
):
    """Construct a configuration dictionary for a single cobra model
    Parameters:
        model_file: string, file path or BiGG Model ID
        name: string, name of the process
        kinetics: dict, kinetic parameters for shared substrates
        reaction_map: dict, maps substrate names to reaction ids
        biomass_identifier: string, name of the biomass reaction
        bounds: dict, bounds for exchange reactions
    Returns:
        config: dict, config dictionary for a single species dFBA
    """
    model = model_from_file(model_file)
    if name is None:
        name = model.id
    if reaction_map is None:
        reaction_map = {
            "glucose": "EX_glc__D_e",
            "acetate": "EX_ac_e"
        }
    if bounds is None:
        bounds = {
            "EX_o2_e": {"lower": -2, "upper": None},
            "ATPM": {"lower": 1, "upper": 1}
        }
    if kinetics is None:
        kinetics = {
            "glucose": (0.02, 15),
            "acetate": (0.5, 7)}
    if biomass_identifier is None:
        biomass_identifier = get_objective_reaction(model=model)

    return {
        "model_file": model_file,
        "name": name,
        "kinetics": kinetics,
        "reaction_map": reaction_map,
        "biomass_identifier": biomass_identifier,
        "bounds": bounds,
    }

def get_single_dfba_spec(
        model_file="textbook",
        name="species",
        config=None,
        interval=1.0
):
    """Constructs a configuration dictionary for a dynamic FBA process
    Parameters:
        model : str, cobra model identifier or path to xml cobra model file
        name: str, identifier for the model, usually species/strain name
        config: dict, config for DFBA Process. If none provided, uses default generated using `dfba_config()`
    Returns:
        dict: dict, specification dictionary for a single species dFBA
    """

    if config is None:
        config = dfba_config(model_file=model_file, name=name)

    return {
        "_type": "process",
        "address": "local:DFBA",
        "config": config,
        "inputs": {
            "shared_environment": ["shared environment"],
            "current_update": ["dFBA Results"]
        },
        "outputs": {
            "dfba_update": ["dFBA Results", name]
        },
        "interval": interval
    }

def environment_spec():
    """Construct spec dictionary for UpdateEnvironment step"""
    return {
        "_type": "process",
        "address": "local:UpdateEnvironment",
        "config": {},
        "inputs": {
            "species_updates": ["dFBA Results"],
            "shared_environment": ["shared environment"]
        },
        "outputs": {
            "shared_environment": ["shared environment"],
        }
    }

def get_chemo_spec(config=None):
    """Constructs a configuration dictionary for the Chemostat process.
    Parameters:
        config: dict, Chemostat configuration dictionary
    Returns:
        dict, spec for Chemostat process
    """
    if config is None:
        raise ValueError("Error: Please provide config")
    return {
        "_type": "process",
        "address": "local:Chemostat",
        "config": config,
        "inputs": {
            "shared_environment": ["shared environment"],
            "global_time": ["global_time"],
        },
        "outputs": {
            "shared_environment": ["shared environment"],
        },
    }

def get_wave_spec(config=None):
    """Constructs a configuration dictionary for the WaveFunction process.
    Parameters:
        config: dict, WaveFunction configuration dictionary
    Returns
        dict, spec for WaveFunction process
    """
    if config is None:
        raise ValueError("Error: Please provide config")
    return {
        "_type": "process",
        "address": "local:WaveFunction",
        "config": config,
        "inputs": {
            "shared_environment": ["shared environment"],
            "global_time": ["global_time"],
        },
        "outputs": {
            "shared_environment": ["shared environment"],
        },
    }

def get_injector_spec(config=None):
    """Constructs a configuration dictionary for the Injector process.
    Parameters:
        config: dict, Injector configuration dictionary
    Returns:
        dict, spec for Injector process
    """
    if config is None:
        raise ValueError("Error: Please provide config")
    return {
        "_type": "process",
        "address": "local:Injector",
        "config": config,
        "inputs": {
            "shared_environment": ["shared environment"],
            "global_time": ["global_time"],
        },
        "outputs": {
            "shared_environment": ["shared environment"],
        },
    }

def run_single_dfba_spec(model_file="textbook"):
    model = model_from_file(model_file)
    exchanges = get_exchanges(model_file=model_file, medium_type='exchange')
    substrates = get_substrates(model_file=model_file, exchanges=exchanges)
    reaction_map = get_reaction_map(model_file=model_file, exchanges=exchanges)
    kinetics = get_kinetics(model_file=model_file, exchanges=exchanges)
    bounds=None
    biomass_identifier = get_objective_reaction(model_file=model_file)
    config = dfba_config(
        model_file=model_file,
        name='E.coli Core',
        kinetics=kinetics,
        reaction_map=reaction_map,
        biomass_identifier=biomass_identifier,
        bounds=bounds,
    )
    spec = get_single_dfba_spec(
        model_file=model_file,
        name='E.coli Core',
        config=config,
        interval=1.0
    )

    pprint.pprint(spec)
    print("")
    print(f"Exchanges: {exchanges}")
    print(f"Substrates: {substrates}")
    print(f"Reaction Map: {reaction_map}")
    print(f"Kinetics: {kinetics}")
    print(f"Biomass Identifier: {biomass_identifier}")
    print(f"Bounds: {bounds}")

def run_initial_counts(model_file="textbook"):
    exchanges = get_exchanges(model_file=model_file, medium_type='exchange')
    initial_counts = get_initial_counts(
        model_dict={'E.coli Core': 'textbook'},
        biomass=0.1,
        initial_value=20,
        exchanges=exchanges
    )

    pprint.pprint(initial_counts)

if __name__ == '__main__':
    # run_single_dfba_spec(model_file="textbook")
    run_initial_counts(model_file="textbook")