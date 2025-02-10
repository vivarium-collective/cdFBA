"""This module contains some methods to obtain the reaction map, initial conditions, and kinetic parameters needed
for dFBA simulations from the minimal medium requirements of the wild-type species.

CAUTION: The initial conditions, and kinetics dataframes provide default parameter values and need to be changed as needed
"""
from copy import deepcopy

#TODO: make notebook demo for functions and make them pure functions

from cobra.io import load_model, read_sbml_model
from cobra.medium import minimal_medium
import pprint
import re
import copy


def make_cdfba_composite(model_dict, kinetic_params=None):
    pass

class DFBAconfig:
 
    def __init__(self, model, medium_type='default'):
        """Creates a medium object.ß

        Parameters
        ----------
        model : cobra model representing the wild-type cell
        medium_type : string, 
            'default' uses the default cobra model medium
            'minimal' uses the minimal medium for the model
            'exchange' uses all exchange fluxes for the model 

        Note
        ----
        Instances of this class provides the reaction mapping, initial conditions
        and kinetic parameters for a COBRA model based on both the minimal and default
        medium of the model. The class methods generate these data within the __init__
        method.
        """
        self.model = model

        if medium_type=='default':
            self.medium = self.model.medium
        if medium_type=='minimal':
            self.medium = minimal_medium(self.model, self.model.slim_optimize()).to_dict()
        if medium_type=='exchange':
            medium = {reaction.id:reaction.upper_bound for reaction in self.model.exchanges}
            medium.update(self.model.medium)
            self.medium = medium
        
        self.substrates = self.get_substrates()
        self.reaction_map = self.get_reaction_map()
        self.kinetics = self.get_kinetics()
        self.biomass_indentifier = get_objective_reaction(self.model)
        
    def get_substrates(self):
        """Returns a list of substrates from the model.
    
        Parameters
        ----------
        medium : DFBAconfig.medium or DFBAconfig.min_medium
        
        Returns
        -------
        substrates : list, list of names of substrates required by the model organism
        """
    
        #obtain substrate names
        substrates = []
        for item in [getattr(self.model.reactions, i).name for i in self.medium.keys()]:
            match = re.match(r"(.*) exchange|exchange reaction for (.*)|Exchange of (.*)", item, re.IGNORECASE)
            if match:
                substrates.append(match.group(1) or match.group(2) or match.group(3))
            else:
                substrates.append(item)
        return substrates
        
    def get_reaction_map(self):
        """Returns a reaction_name_map dictionary from a medium dictionary as obtained
        from model.medium or cobra.medium.minimum_medium()
        
        Parameters
        ----------
        medium : DFBAconfig.medium or DFBAconfig.min_medium
        substrates : list, list of names of substrates required by the model organism
        
        Returns
        -------
        reaction_name_map : dict, maps substrate names to reactions
        """    
        substrates = copy.deepcopy(self.substrates)
        ids = copy.deepcopy(list(self.medium.keys()))
        reaction_name_map = {}
        for i in range(len(substrates)):
            reaction_name_map[substrates[i]] = ids[i]
        return reaction_name_map
    
    def get_kinetics(self):
        """Returns default kinetic parameters dictionary
        Values are tuples of the form (km, vmax)
        
        Parameters
        ----------
        substrates    : list, list of names of substrates required by the model organism
        """  
        kinetics = {key: (0.5, 2.0) for key in self.substrates}
        
        return kinetics


def get_objective_reaction(model):
    """get a string with the name of the objective function of a cobra model

    Parameters:
    -----------
    model: cobrapy model

    Returns:
    --------
    objective_reaction: string, name of the objective reaction (biomass reaction by default)
    """

    expression = f"{model.objective.expression}"
    match = re.search(r'1\.0\*([^\s]+)', expression)

    if match:
        objective_reaction = match.group(1)

    return objective_reaction

def initial_conditions(model, biomass=0.1, factor=1.0, medium_type='default', name=None, default_concentration = None):
    """Returns an initial condition dict based on medium
    
    Parameters
    ----------
    model: string, cobrapy model name
    substrates : list, list of names of substrates required by the model organism
    biomass : float, initial biomass for all species
    factor : float, factor to multiply minimum medium concentrations
    
    Returns
    -------
    conditions : dict, initial conditions dictionary
    """ 
    medium = DFBAconfig(model=model, medium_type=medium_type)
    conditions = {}
    if name is None:
        conditions.update({get_objective_reaction(model): biomass})
    else:
        conditions.update({name: biomass})
    
    substrates = medium.substrates

    if default_concentration is not None:
        substrate_values = dict(zip(substrates, list(medium.medium.values())))
    else:
        substrate_values = dict(zip(substrates, list(medium.medium.values())))


    for key in substrate_values:
        substrate_values[key] *= factor
    
    conditions.update(substrate_values)
    
    return conditions

def model_from_file(model_file='textbook'):
    if not "xml" in model_file:
            # use the textbook model if no model file is provided
            # TODO: Also handle JSON or .mat model files
            model = load_model(model_file)
    elif isinstance(model_file, str):
        model = read_sbml_model(model_file)
    else:
        # error handling
        raise ValueError("Invalid model file")
    return model

def model_list(model_files=None):
    """generates list of cobra models from files/model ids.

    Parameters:
    -----------
    model_files: list, list of strings with the model file paths or model ids

    Returns:
    --------
    model_list: list, list of cobra models
    """

    return [model_from_file(model_file) for model_file in model_files]

def dfba_config(
        model_file="textbook",
        name=None,
        kinetics=None,
        reaction_map=None,
        biomass_identifier=None,
        bounds=None
):
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
        "interval": 0.1,
    }

def dfba_config_from_model(
        model_file="textbook",
        medium_type ="default",
        name=None,
        bounds=None
):
    model = model_from_file(model_file=model_file)
    if name is None:
        name = model.id
    if bounds is None:
        bounds = {
            "EX_o2_e": {"lower": -2, "upper": None},
            "ATPM": {"lower": 1, "upper": 1}
        }

    dfbaconfig = DFBAconfig(model, medium_type=medium_type)
    kinetics = dfbaconfig.kinetics
    reaction_map = dfbaconfig.reaction_map
    biomass_identifier = dfbaconfig.get_objective_reaction(model)

    return {
        "model": model,
        "name": name,
        "kinetics": kinetics,
        "reaction_map": reaction_map,
        "biomass_identifier": biomass_identifier,
        "bounds": bounds
    }

def get_single_dfba_spec(
        model_file="textbook",
        name="species",
        config=None #TODO: add kwargs
):
    """
    Constructs a configuration dictionary for a dynamic FBA process with optional path indices.

    This function builds a process specification for use with a dynamic FBA system. It allows
    specification of substrate molecule IDs and optionally appends indices to the paths for those substrates.

    Parameters:
    -----------
    model : str, cobra model identifier or path to xml cobra model file
    name: str, identifier for the model, usually species/strain name
    config: dict, config for DFBA Process. If none provided, uses default generated using `dfba_config()`

    Returns:
    --------
    dict: A dictionary containing the process type, address, configuration, and paths for inputs
        and outputs based on the specified molecule IDs and indices.
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
        "interval": 1.0
    }

def environment_spec():
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

def community_dfba_spec(
        species_list = None,
        from_model=False,
        medium_type='default'
):
    stores = {
        'shared environment': 'any',
        'dFBA Results': 'any',
    }

    dfba_processes = {}

    if from_model:
        for model in species_list:
            dfba_processes.update(

            )
    #TODO: Finish this function

def initial_environment(volume=1, initial_counts=None, species_list=None):
    """
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

def get_chemo_spec(config=None):
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

def run_config(medium_type='default'):
    # model = model_from_file('/Users/tasnifrahman/Research/ecolicommunity/agora2/Eubacterium_rectale_ERR1203958.xml')
    model = load_model('iAF1260')
    config = DFBAconfig(model, medium_type=medium_type)
    # pprint.pprint(config.substrates)
    pprint.pprint(config.reaction_map)

if __name__ == '__main__':
    run_config(medium_type='exchange')
