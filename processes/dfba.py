import cobra
from cobra.io import load_model
from cobra.medium import minimal_medium
from process_bigraph import Process, Step, Composite

model = load_model("iAF1260")


class DFBA(Process):
    """Performs single time-step of dynamic FBA
    
    Parameters:
    -----------
    model_file: string, math to cobra model file  
    """
    config_schema = {
        'model_file': 'string',
        'kinetics': 'map[tuple[float, float]]',
        'reaction_map': 'map[string]',
        'biomass_identifier': 'string',
        'bounds': 'map[bounds]',
    }
    
    def __init__(self, config, core):
        super().__init__(config, core)
        
        if not "xml" in self.config["model_file"]:
            # use the textbook model if no model file is provided
            # TODO: Also handle JSON or .mat model files
            self.model = load_model(self.config["model_file"])
        elif isinstance(self.config["model_file"], str):
            self.model = cobra.io.read_sbml_model(self.config["model_file"])
        else:
            # error handling
            raise ValueError("Invalid model file")
        if self.config['bounds'] is not None:
            for reaction_id, bounds in self.config["bounds"].items():
                if bounds["lower"] is not None:
                    self.model.reactions.get_by_id(reaction_id).lower_bound = bounds["lower"]
                if bounds["upper"] is not None:
                    self.model.reactions.get_by_id(reaction_id).upper_bound = bounds["upper"]

    def inputs(self):
        return {
             'shared_environment': 'any' #initial conditions for time-step
        }
        
    def outputs(self):
        return {
             'dfba_update': 'any'
        }
        
    def update(self, inputs, interval):
            
        current_state = inputs['shared_environment'].copy()
        state_update = inputs['shared_environment'].copy()
    
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
        state_update[self.config["model_file"]] = biomass_growth_rate * current_biomass * interval
    
        ## update substrates
        for substrate_id, reaction_id in self.config["reaction_map"].items():
            flux = solution.fluxes[reaction_id]
            current_substrate_conc = current_state[substrate_id]
            state_update[substrate_id] = flux * current_biomass * interval

        return {'dfba_update': state_update}
    
class UpdateEnvironment(Step):
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)
        
    def inputs(self):
        return {
             'substrates': 'map[positive_float]' #initial conditions for time-step
        }
        
    def outputs(self):
        return {
             'port2': 'any'
        }

    def update(self, inputs, interval):

        species_updates = inputs['species_updates']
        shared_environment = inputs['shared_environment']
        update = {}
        return {'port2': update}

def dfba_config(
        model_file="textbook",
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
        "model_file": model_file,
        "kinetics": kinetics,
        "reaction_map": reaction_map,
        "biomass_identifier": biomass_identifier,
        "bounds": bounds
    }

def dfba_config_from_model(
        model_file="textbook",
        kinetics=None,
        reaction_map=None,
        biomass_identifier="biomass",
        bounds=None
):
#TODO: finish this method
    return {
        "model_file": model_file,
        "kinetics": kinetics,
        "reaction_map": reaction_map,
        "biomass_identifier": biomass_identifier,
        "bounds": bounds
    }
    
def get_single_dfba_spec(
        model_file="textbook",
        name="species",
        config=None
):
    """
    Constructs a configuration dictionary for a dynamic FBA process with optional path indices.

    This function builds a process specification for use with a dynamic FBA system. It allows
    specification of substrate molecule IDs and optionally appends indices to the paths for those substrates.

    Parameters:
        mol_ids (list of str, optional): List of molecule IDs to include in the process. Defaults to
                                         ["glucose", "acetate", "biomass"].
        path (list of str, optional): The base path to prepend to each molecule ID. Defaults to ["..", "fields"].
        i (int, optional): The first index to append to the path for each molecule, if not None.
        j (int, optional): The second index to append to the path for each molecule, if not None.

    Returns:
        dict: A dictionary containing the process type, address, configuration, and paths for inputs
              and outputs based on the specified molecule IDs and indices.
    """

    if config is None:
        config = dfba_config(model_file=model_file)

    return {
        "_type": "process",
        "address": "local:DynamicFBA",
        "config": dfba_config(model_file=model_file),
        "inputs": {
            "read environment": "shared environment"
        },
        "outputs": {
            "substrate updates": ["dFBA Results", f"{name}"]
        }
    }

#TODO write config and composite spec generator functions
#TODO write update_environment Step