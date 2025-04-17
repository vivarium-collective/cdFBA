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

#global variables
SHARED_ENVIRONMENT = "Shared Environment"
SPECIES_STORE = "Species"
DFBA_RESULTS = "dFBA Results"
THRESHOLDS = "Thresholds"

#basic functions
def model_from_file(model_file="textbook"):
    """Returns a cobra model from a model file path or BiGG Model ID
    Parameters:
        model_file: str, file path or BiGG Model ID
    Returns:
        model: cobra model
    """
    #check for model type and load model
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

def get_model_dict(model_dict):
    """
    Returns a dictionary of the cobra model from a dictionary of model IDs/file paths
    Args:
        model_dict: dict, dictionary of BIGG IDs or model paths
    Returns:
        model_dict: dict, dictionary of cobra models
    """
    return {key:model_from_file(value) for key, value in model_dict.items()}

#set value functions
def set_counts(spec, counts):
    """Set counts for given metabolites in the shared environment
        Parameters:
            spec : dict, cdFBA specification dictionary
            counts : dict, dictionary with substrate names as keys and counts as values
        """
    for substrate, count in counts.items():
        #make sure substrates are present in the cdfba model
        if substrate not in spec[SHARED_ENVIRONMENT]["concentrations"].keys():
            raise ValueError(f"{substrate} is not in shared environment")
        else:
            spec[SHARED_ENVIRONMENT]["counts"][substrate] = count
            spec[SHARED_ENVIRONMENT]["concentrations"][substrate] = spec[SHARED_ENVIRONMENT]["counts"][substrate]/spec[SHARED_ENVIRONMENT]["volume"]

def set_concentration(spec, concentrations):
    """Set concentration for given metabolites in the shared environment
    Parameters:
        spec : dict, cdFBA specification dictionary
        concentrations : dict, dictionary with substrate names as keys and concentrations as values
    """
    # make sure substrates are present in the cdfba model
    for substrate, concentration in concentrations.items():
        if substrate not in spec[SHARED_ENVIRONMENT]["concentrations"].keys():
            raise ValueError(f"{substrate} is not in shared environment")
        else:
            spec[SHARED_ENVIRONMENT]["concentrations"][substrate] = concentration
            spec[SHARED_ENVIRONMENT]["counts"][substrate] = spec[SHARED_ENVIRONMENT]["concentrations"][substrate] * spec[SHARED_ENVIRONMENT]["volume"]

def set_kinetics(species, spec, kinetics):
    """Set bounds for given metabolites and species
    Parameters:
        species: str, name of species
        spec : dict, cdFBA specification dictionary
        kinetics : dict, substrate names as keys and kinetics parameters in tuples as values (Km, Vmax)
    """
    for substrate, kinetic_params in kinetics.items():
        if substrate not in spec[SHARED_ENVIRONMENT]["concentrations"].keys():
            raise ValueError(f"{substrate} is not in shared environment")
        else:
            spec[SPECIES_STORE][species]["config"]["kinetics"][substrate] = kinetic_params

#single species functions
def get_exchanges(model_file="textbook", medium_type="exchange"):
    """
    Parameters:
        model_file: str, file path or BiGG Model ID, OR
                    cobra model
        medium_type:
            "default" uses the default cobra model medium
            "minimal" uses the minimal medium for the model
            "exchange" uses all exchange fluxes for the model
            defaults to "exchange"
    Returns:
        exchanges: list of exchange reaction IDs
    """
    if isinstance(model_file, str):
        model = model_from_file(model_file)
    else:
        model = model_file

    if medium_type == "default":
        medium = model.medium
    if medium_type == "minimal":
        medium = minimal_medium(model, model.slim_optimize()).to_dict()
    if medium_type == "exchange":
        medium = {reaction.id: reaction.upper_bound for reaction in model.exchanges}
        medium.update(model.medium)
        medium = medium
    if not medium_type in ["default", "minimal", "exchange"]:
        raise ValueError("Invalid medium type")

    return list(medium.keys())

def get_substrates(model_file="textbook", exchanges=None):
    """Returns a list of substrates from the model.
    Parameters:
    model_file : str, file path or BiGG Model ID
    exchanges : lst, a list of exchange reaction ids
    Returns:
    substrates : lst, list of names of substrates required by the model organism
    """
    if isinstance(model_file, str):
        model = model_from_file(model_file)
    else:
        model = model_file
    if exchanges is None:
        exchanges = get_exchanges(model) #gets all exchange reactions by default
    substrates = [item for item in [list(getattr(model.reactions, i).metabolites.keys())[0].name for i in exchanges if hasattr(model.reactions, i)]]
    return substrates
        
def get_reaction_map(model_file="textbook", exchanges=None):
    """Returns a reaction_name_map dictionary from a medium dictionary as obtained
    from model.medium or cobra.medium.minimum_medium()
    Parameters:
        model_file : str, file path or BiGG Model ID
        exchanges : lst, list of names of substrates required by the model organism
    Returns:
        reaction_name_map : dict, maps substrate names to reactions
    """
    if isinstance(model_file, str):
        model = model_from_file(model_file)
    else:
        model = model_file
    if exchanges is None:
        exchanges = get_exchanges(model)
    substrates = get_substrates(model, exchanges)
    ids = exchanges
    reaction_name_map = {}
    for i in range(len(substrates)):
        reaction_name_map[substrates[i]] = ids[i]
    return reaction_name_map
    
def get_kinetics(model_file="textbook", exchanges=None):
    """Returns default kinetic parameters dictionary. Values are tuples of the form (km, vmax)
    Parameters:
        model_file : str, file path or BiGG Model ID
        exchanges : lst, list of exchange reaction ids
    """
    if isinstance(model_file, str):
        model = model_from_file(model_file)
    else:
        model = model_file
    if exchanges is None:
        exchanges = get_exchanges(model)
    kinetics = {key: (0.5, 2.0) for key in get_substrates(model, exchanges)}
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
    return {reaction_map[key]: {"lower": lower, "upper": upper} for key in reaction_map.keys()}

def get_objective_reaction(model_file = "textbook"):
    """get a string with the name of the objective function of a cobra model
    Parameters:
        model_file: str, BiGG Model ID or model path OR
                    cobra model instance
    Returns:
        objective_reaction: str, name of the objective reaction (biomass reaction by default)
    """
    if isinstance(model_file, str):
        model = model_from_file(model_file)
    else:
        model = model_file
    expression = f"{model.objective.expression}"
    match = re.search(r"1\.0\*([^\s]+)", expression)

    if match:
        objective_reaction = match.group(1)

    return objective_reaction

def dfba_config(
        model_file="textbook",
        model = None,
        name=None,
        kinetics=None,
        reaction_map=None,
        bounds=None,
        changes=None,
):
    """Construct a configuration dictionary for a single cobra model
    Parameters:
        model_file: str, file path or BiGG Model ID
        model: cobra model (optional, provide if model already loaded to avoid reloading)
        name: str, name of the process
        kinetics: dict, kinetic parameters for shared substrates
        reaction_map: dict, maps substrate names to reaction ids
        bounds: dict, bounds for exchange reactions
        changes: dict, changes to apply to the model
    Returns:
        config: dict, config dictionary for a single species dFBA
    """
    if model is None:
        model = model_from_file(model_file)
    if name is None:
        name = model.id
    if reaction_map is None:
        reaction_map = {
            "D-Glucose": "EX_glc__D_e",
            "Acetate": "EX_ac_e"
        }
    if bounds is None:
        bounds = {
            "EX_o2_e": {"lower": -2, "upper": None},
            "ATPM": {"lower": 1, "upper": 1}
        }
    if kinetics is None:
        kinetics = {
            "glucose": (0.02, 15),
            "acetate": (0.5, 7)
        }
    if changes is None:
        changes = {
            "gene_knockout": [],
            "reaction_knockout": [],
            "bounds": {},
            "kinetics": {},
        }
    return {
        "model_file": model_file,
        "name": name,
        "kinetics": kinetics,
        "reaction_map": reaction_map,
        "bounds": bounds,
        "changes": changes,
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
        config: dict, config for dFBA Process. If none provided, uses default generated using `dfba_config()`
        interval: float, interval between consecutive dFBA calculations
    Returns:
        dict: dict, specification dictionary for a single species dFBA
    """
    if isinstance(model_file, str):
        model = model_from_file(model_file)
    else:
        model = model_file
    if config is None:
        config = dfba_config(model_file=model, name=name)

    return {
        "_type": "process",
        "address": "local:dFBA",
        "config": config,
        "inputs": {
            "shared_environment": ["..", SHARED_ENVIRONMENT],
            "current_update": ["..", DFBA_RESULTS]
        },
        "outputs": {
            "dfba_update": ["..", DFBA_RESULTS, name]
        },
        "interval": interval
    }

#multi-species functions
def make_cdfba_composite(model_dict, medium_type=None, exchanges=None, volume=1, interval=1.0):
    """Construct a cdfba composite spec with all exhange metabolites included.
    Parameters:
        model_dict : dict, dictionary with cdfba process names as keys and model name/path as values
        medium_type : str/lst, if str, pick one of:
            "default" uses the default cobra model medium
            "minimal" uses the minimal medium for the model
            "exchange" uses all exchange fluxes for the model
            MUST be None if exchanges is provided
        exchanges: a list of exchange reaction ids. MUST be None if medium_type is provided
        volume: float, volume of cdfba composite
        interval: float, interval between consecutive dFBA calculations
    Returns:
        spec : dict, cdfba composite spec
    """
    #load models
    models_dict = get_model_dict(model_dict)
    #initialize spec
    spec = {DFBA_RESULTS: {}}
    #ensure only one of medium_type or exchanges is provided
    if medium_type is None:
        if exchanges is None:
            raise ValueError("Must provide medium_type or exchanges list")
    if medium_type is not None:
        if exchanges is not None:
            raise ValueError("Provide only on of medium_type or exchanges list")
    #get union of exchange reactions from all species if exchanges not provided
    if exchanges is None:
        env_exchanges = get_combined_exchanges(models_dict, medium_type=medium_type)
    else:
        env_exchanges = exchanges
    #set initial environment
    initial_counts = get_initial_counts(models_dict, exchanges=env_exchanges)
    initial_env = initial_environment(volume=volume, initial_counts=initial_counts, species_list=models_dict.keys())
    spec[SHARED_ENVIRONMENT] = initial_env
    #generate all dFBA processes
    spec[SPECIES_STORE] = {}
    for model_name, model_file in models_dict.items():
        #get model-specific exchanges if exchanges not provided
        if exchanges is None:
            model_exchanges = get_exchanges(model_file=model_file, medium_type=medium_type)
            model_exchanges = [exchange for exchange in model_exchanges if exchange in env_exchanges]
        else:
            model_exchanges = env_exchanges
        #get list of substrates
        substrates = get_substrates(model_file=model_file, exchanges=model_exchanges)
        #get default kinetics parameters
        kinetics = get_kinetics(model_file=model_file, exchanges=model_exchanges)
        #get reaction map
        reaction_map = get_reaction_map(model_file=model_file, exchanges=model_exchanges)
        #set default bounds
        bounds = {}

        config = dfba_config(
            model_file=model_dict[model_name],
            name=model_name,
            kinetics=kinetics,
            reaction_map=reaction_map,
            bounds=bounds
        )
        model_spec = get_single_dfba_spec(
            model_file=model_file,
            name=model_name,
            config=config,
            interval=interval
        )
        #add dFBA spec to composite spec
        spec[SPECIES_STORE][model_name] = model_spec
        #initialize dFBA results store
        spec[DFBA_RESULTS][model_name] = {substrate: 0 for substrate in substrates}
        spec[DFBA_RESULTS][model_name].update({model_name: 0})
    #add UpdateEnvironment step spec
    spec["update environment"] = environment_spec()
    return spec

def get_combined_exchanges(model_dict, medium_type=None):
    """Returns a list of exchange reaction ids for multiple species - containing every
    Parameters:
        model_dict: dict, dictionary with cdfba process names as keys and model name/path as values
        medium_type: str/lst,pick one of:
        "default" uses the default cobra model medium
        "minimal" uses the minimal medium for the model
        "exchange" uses all exchange fluxes for the model
    Returns:
        env_exchanges: list, list of
    """
    env_exchanges = []
    for name, model_file in model_dict.items():
        env_exchanges.extend(get_exchanges(model_file=model_file, medium_type=medium_type))
    env_exchanges = list(set(env_exchanges))
    return env_exchanges

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
    if exchanges is None:
        raise ValueError("Must provide list of exchange reaction ids")
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
            "species_updates": [DFBA_RESULTS],
            "shared_environment": [SHARED_ENVIRONMENT]
        },
        "outputs": {
            "counts": [SHARED_ENVIRONMENT, "counts"],
        }
    }

def get_static_spec(config=None, interval=1.0):
    """Constructs a configuration dictionary for the StaticConcentration process.
    Parameters:
        config: dict, StaticConcentration configuration dictionary
    Returns:
        dict, spec for StaticConcentration process
    """
    if config is None:
        raise ValueError("Error: Please provide config")
    return {
        "_type": "process",
        "address": "local:StaticConcentration",
        "config": config,
        "inputs": {
            "shared_environment": [SHARED_ENVIRONMENT],
            "global_time": ["global_time"],
        },
        "outputs": {
            "shared_environment": [SHARED_ENVIRONMENT],
        },
        "interval": interval,
    }

def get_wave_spec(config=None, interval=1.0):
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
            "shared_environment": [SHARED_ENVIRONMENT],
            "global_time": ["global_time"],
        },
        "outputs": {
            "shared_environment": [SHARED_ENVIRONMENT],
        },
        "interval": interval,
    }

def get_injector_spec(config=None, interval=1.0):
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
            "shared_environment": [SHARED_ENVIRONMENT],
            "global_time": ["global_time"],
        },
        "outputs": {
            "shared_environment": [SHARED_ENVIRONMENT],
        },
        interval: interval,
    }

#=======
# TESTS
#=======

def run_single_dfba_spec(model_file="textbook"):
    model = model_from_file(model_file)
    exchanges = get_exchanges(model_file=model_file, medium_type="exchange")
    substrates = get_substrates(model_file=model_file, exchanges=exchanges)
    reaction_map = get_reaction_map(model_file=model_file, exchanges=exchanges)
    kinetics = get_kinetics(model_file=model_file, exchanges=exchanges)
    bounds=None
    biomass_identifier = get_objective_reaction(model_file=model_file)
    config = dfba_config(
        model_file=model_file,
        name="E.coli Core",
        kinetics=kinetics,
        reaction_map=reaction_map,
        biomass_identifier=biomass_identifier,
        bounds=bounds,
    )
    spec = get_single_dfba_spec(
        model_file=model_file,
        name="E.coli Core",
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
    exchanges = get_exchanges(model_file=model_file, medium_type="exchange")
    initial_counts = get_initial_counts(
        model_dict={"E.coli Core": "textbook"},
        biomass=0.1,
        initial_value=20,
        exchanges=exchanges
    )

    pprint.pprint(initial_counts)

def run_cdfba_spec():
    model_dict = {
        "E.coli" : "iAF1260",
        "S.flexneri" : "iSFxv_1172"
    }

    pprint.pprint(make_cdfba_composite(
        model_dict=model_dict,
        medium_type="minimal",
        exchanges=None,
        volume=1,
        interval=1.0
        )
    )

if __name__ == "__main__":
    # run_single_dfba_spec(model_file="textbook")
    # run_initial_counts(model_file="textbook")
    run_cdfba_spec()
