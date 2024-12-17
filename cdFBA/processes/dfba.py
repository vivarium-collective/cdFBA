import sys
import os 
import cobra
import random
from cobra.io import load_model
from cobra.medium import minimal_medium
from process_bigraph import Process, Step, Composite

sys.path.append('/Users/tasnifrahman/Research/cdFBA')

from cdFBA.utils import DFBAconfig


class DFBA(Process):
    """Performs single time-step of dynamic FBA
    
    Parameters:
    -----------
    model_file: string, math to cobra model file  
    """
    config_schema = {
        "model_file": "any",
        "kinetics": "any",
        "reaction_map": "any",
        "biomass_identifier": "any",
        "bounds": "any",
    }
    
    def __init__(self, config, core):
        super().__init__(config, core)



        if self.config["bounds"] is not None:
            for reaction_id, bounds in self.config["bounds"].items():
                if bounds["lower"] is not None:
                    self.model.reactions.get_by_id(reaction_id).lower_bound = bounds["lower"]
                if bounds["upper"] is not None:
                    self.model.reactions.get_by_id(reaction_id).upper_bound = bounds["upper"]

    def inputs(self):
        return {
             "shared_environment": "any" #initial conditions for time-step
        }
        
    def outputs(self):
        return {
             "dfba_update": "any"
        }
        
    def update(self, inputs, interval):
            
        current_state = inputs["shared_environment"].copy()
        state_update = inputs["shared_environment"].copy()
    
        for substrate_id, reaction_id in self.config["reaction_map"].items():
            Km, Vmax = self.config["kinetics"][substrate_id]
            substrate_concentration = current_state[substrate_id]
            
            # calculate michaelis-menten flux
            flux = Vmax * substrate_concentration / (Km + substrate_concentration)
    
            # use the flux to constrain fba
            self.model.reactions.get_by_id(reaction_id).lower_bound = -flux
    
        # solve fba under these constraints
        solution = self.model.optimize()
    
        # gather the results
        ## update biomass
        biomass_growth_rate = solution.fluxes[self.config["biomass_identifier"]]
        current_biomass = current_state[self.config["biomass_identifier"]]
        state_update[self.config["model"]] = biomass_growth_rate * current_biomass * interval
    
        ## update substrates
        for substrate_id, reaction_id in self.config["reaction_map"].items():
            flux = solution.fluxes[reaction_id]
            current_substrate_conc = current_state[substrate_id]
            state_update[substrate_id] = flux * current_biomass * interval

        return {"dfba_update": state_update}
    
def dfba_config(
        model=load_model("textbook"),
        kinetics=None,
        reaction_map=None,
        biomass_identifier="biomass",
        bounds=None
):
    if reaction_map is None:
        reaction_map = {
            "glucose": "EX_glc__D_e",
            "acetate": "EX_ac_e"}
    if bounds is None:
        bounds = {}
    if kinetics is None:
        kinetics = {
            "glucose": (0.5, 1),
            "acetate": (0.5, 2)}
    return {
        "model": model,
        "kinetics": kinetics,
        "reaction_map": reaction_map,
        "biomass_identifier": biomass_identifier,
        "bounds": bounds
    }

def dfba_config_from_model(
        model=load_model("textbook"),
        medium_type ="default",
        bounds={}
):
    dfbaconfig = DFBAconfig(model, medium_type=medium_type)
    kinetics = dfbaconfig.kinetics
    reaction_map = dfbaconfig.reaction_map
    biomass_identifier = dfbaconfig.get_objective_reaction(model)

    return {
        "model": model,
        "kinetics": kinetics,
        "reaction_map": reaction_map,
        "biomass_identifier": biomass_identifier,
        "bounds": bounds
    }
    
def get_single_dfba_spec(
        model=load_model("textbook"),
        name="species",
        config=None
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
    dict: A dictionary containing the process type, address, configuration, and paths for inputs
        and outputs based on the specified molecule IDs and indices.
    """

    if config is None:
        config = dfba_config(model=model)

    return {
        "_type": "process",
        "address": "local:DFBA",
        "config": config,
        "inputs": {
            "shared_environment": "shared environment"
        },
        "outputs": {
            "dfba_update": ["dFBA Results", f"{name}"]
        }
    }

class UpdateEnvironment(Step):
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)
        
    def inputs(self):
        return {
             "shared_environment": "any",
             "species_updates": "any"
        }
        
    def outputs(self):
        return {
             "shared_environment": "any"
        }

    def update(self, inputs):

        species_updates = inputs["species_updates"]
        shared_environment = inputs["shared_environment"]

        species_list = random.shuffle([species for species in species_updates])
        update = shared_environment.copy()

        for species in species_list:
            update = {key:update[key] + species_updates[species][key] for key in update}

        update = {}
        return {"shared_environment": update}
    
def environment_spec():
    return {
        "_type": "step",
        "address": "local:UpdateEnvironment",
        "config": {},
        "inputs": {
            "species_updates": "dFBA Results",
            "shared_environment": "shared environment"
        },
        "outputs": {
            "shared_environment": "shared environment"
        }
    }

def community_dfba_spec(species_list = [], from_model=False, medium_type='default'):
    stores = {
        'shared environment': 'any',
        'dFBA Results': 'any',
    }

    dfba_processes = {}
    
    if from_model:
        for model in species_list:
            dfba_processes.update(
                
            )

if __name__ == "__main__":
    print(get_single_dfba_spec())