"""This module contains some methods to obtain the reaction map, initial conditions, and kinetic parameters needed
for dFBA simulations from the minimal medium requirements of the wild-type species.

CAUTION: The initial conditions, and kinetics dataframes provide default parameter values and must be changed as needed
CAUTION: Substrate names are different in BiGG and AGORA databases. These functions will not work with two models form
         different sources
"""
from cobra.io import load_model, read_sbml_model, load_json_model, load_yaml_model, load_matlab_model
from cobra.medium import minimal_medium
import pprint
import re

#single species functions
def model_from_file(model_file='textbook'):
    """Returns a cobra model from a model file path or BiGG Model ID
    Parameters:
        model_file: str, file path or BiGG Model ID
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
    """

    Parameters:
        model_file: str, file path or BiGG Model ID
        medium_type:
            'default' uses the default cobra model medium
            'minimal' uses the minimal medium for the model
            'exchange' uses all exchange fluxes for the model
            defaults to 'exchange'
    Returns:
        exchanges: list of exchange reaction IDs
    """
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
    model_file : str, file path or BiGG Model ID
    exchanges : lst, a list of exchange reaction ids
    Returns:
    substrates : lst, list of names of substrates required by the model organism
    """
    if exchanges is None:
        exchanges = get_exchanges(model_file)
    model = model_from_file(model_file)
    substrates = []

    for item in [getattr(model.reactions, i).name for i in exchanges if hasattr(model.reactions, i)]:

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
        model_file : str, file path or BiGG Model ID
        exchanges : lst, list of names of substrates required by the model organism
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
    """Returns default kinetic parameters dictionary. Values are tuples of the form (km, vmax)
    Parameters:
        model_file : str, file path or BiGG Model ID
        exchanges : lst, list of exchange reaction ids
    """
    if exchanges is None:
        exchanges = get_exchanges(model_file)
    model = model_from_file(model_file)
    kinetics = {key: (0.5, 2.0) for key in get_substrates(model_file, exchanges)}
    return kinetics

def get_bounds(reaction_map, upper=1000, lower=-1000):
    """Return dict of default upper and lower bounds for each substrate exchange reaction
    Parameters:
        reaction_map : dict, maps substrate names to reactions
        upper : float, upper bound, default 1000
        lower : float, lower bound, default -1000
    Returns:
        bounds : dict, default bounds for all exchange reactions
    """
    return {reaction_map[key]: {'lower': lower, 'upper': upper} for key in reaction_map.keys()}

def get_objective_reaction(model_file = 'textbook'):
    """get a string with the name of the objective function of a cobra model
    Parameters:
        model: cobrapy model
    Returns:
        objective_reaction: str, name of the objective reaction (biomass reaction by default)
    """
    model = model_from_file(model_file)
    expression = f"{model.objective.expression}"
    match = re.search(r'1\.0\*([^\s]+)', expression)

    if match:
        objective_reaction = match.group(1)

    return objective_reaction

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
        model_file: str, file path or BiGG Model ID
        name: str, name of the process
        kinetics: dict, kinetic parameters for shared substrates
        reaction_map: dict, maps substrate names to reaction ids
        biomass_identifier: str, name of the biomass reaction
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
        biomass_identifier = get_objective_reaction(model_file=model_file)

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
        model_file: str, file path or BiGG Model ID
        name: str, identifier for the model, usually species/strain name
        config: dict, config for DFBA Process. If none provided, uses default generated using `dfba_config()`
        interval: float, interval between consecutive dFBA calculations
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

#multi-species functions
def make_cdfba_composite(model_dict, medium_type=None, exchanges=None, volume=1, interval=1.0):
    """Construct a cdfba composite spec with all exhange metabolites included.
    Parameters:
        model_dict : dict, dictionary with cdfba process names as keys and model name/path as values
        medium_type : str/lst, if str, pick one of:
            'default' uses the default cobra model medium
            'minimal' uses the minimal medium for the model
            'exchange' uses all exchange fluxes for the model
            MUST be None if exchanges is provided
        exchanges: a list of exchange reaction ids. MUST be None if medium_type is provided
        volume: float, volume of cdfba composite
        interval: float, interval between consecutive dFBA calculations
    Returns:
        spec : dict, cdfba composite spec
    """
    spec = {'dFBA Results': {}}
    if medium_type is None:
        if exchanges is None:
            raise ValueError("Must provide medium_type or exchanges list")

    if medium_type is not None:
        if exchanges is not None:
            raise ValueError("Provide only on of medium_type or exchanges list")

    if exchanges is None:
        env_exchanges = []
        for name, model_file in model_dict.items():
            env_exchanges.extend(get_exchanges(model_file=model_file, medium_type=medium_type))

        env_exchanges = list(set(env_exchanges))
    else:
        env_exchanges = exchanges

    initial_counts = get_initial_counts(model_dict, exchanges=env_exchanges)
    initial_env = initial_environment(volume=volume, initial_counts=initial_counts, species_list=model_dict.keys())
    spec['shared environment'] = initial_env

    for model_name, model_file in model_dict.items():
        if exchanges is None:
            model_exchanges = get_exchanges(model_file=model_file, medium_type=medium_type)
        else:
            model_exchanges = exchanges
        substrates = get_substrates(model_file=model_file, exchanges=model_exchanges)
        kinetics = get_kinetics(model_file=model_file, exchanges=model_exchanges)
        reaction_map = get_reaction_map(model_file=model_file, exchanges=model_exchanges)
        biomass_identifier = get_objective_reaction(model_file=model_file)
        bounds = {}

        config = dfba_config(
            model_file=model_file,
            name=model_name,
            kinetics=kinetics,
            reaction_map=reaction_map,
            biomass_identifier=biomass_identifier,
            bounds=bounds
        )
        model_spec = get_single_dfba_spec(model_file=model_file, name=model_name, config=config, interval=interval)
        spec[model_name] = model_spec

        spec['dFBA Results'][model_name] = {substrate: 0 for substrate in substrates}
        spec['dFBA Results'][model_name].update({model_name: 0})
    spec['update environment'] = environment_spec()
    return spec

def get_initial_counts(model_dict, biomass=0.5, initial_value=20, exchanges=None):
    """Returns an initial condition dict based on medium
    Parameters:
        model_dict: dict, dictionary with cdfba process names as keys and model name/path as values
        biomass : float, initial biomass for all species
        initial_value : float, initial counts of all species
        exchanges: lst, list of exchange reaction ids
    Returns:
        conditions : dict, initial conditions dictionary
    """
    all_substrates = []
    for model_name in model_dict.keys():
        model_file=model_dict[model_name]
        substrates = get_substrates(model_file, exchanges)
        all_substrates.extend(substrates)
    all_substrates = list(set(all_substrates))
    conditions = {substrate:initial_value for substrate in all_substrates}
    biomasses = {model:biomass for model in model_dict.keys()}
    conditions = conditions | biomasses
    return conditions

def initial_environment(volume=1, initial_counts=None, species_list=None):
    """Construct initial shared environment store
    Parameters:
        volume : float, volume of the environment
        initial_counts : dict, initial counts of each substrate and species biomass in the environment
        species_list : list of strings, list of dfba species names)
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

#environmental process/step related functions
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

def run_cdfba_spec():
    model_dict = {
        'E.coli' : 'iAF1260',
        'S. flexneri' : 'iSFxv_1172'
    }

    pprint.pprint(make_cdfba_composite(
        model_dict=model_dict,
        medium_type='exchange',
        exchanges=None,
        volume=1,
        interval=1.0
        )
    )

if __name__ == '__main__':
    # run_single_dfba_spec(model_file="textbook")
    # run_initial_counts(model_file="textbook")
    run_cdfba_spec()