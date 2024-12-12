from cobra.io import load_model
from cobra.medium import minimal_medium
from process_bigraph import Process, ProcessTypes, Composite

model = load_model("iAF1260")


class dFBA(Process):
    """Performs single time-step of dynamic FBA
    
    Parameters:
    -----------
    model: 
    """
    config_schema = {
        'model_file': 'string',
        'kinetics': 'map[tuple[float, float]]',
        'reaction_map': 'map[string]',
        'biomass_identifier': 'list',
    }
    
    def __init__(self, config, core):
        super().__init__(config, core)
        
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
            model.reactions.get_by_id(reaction_id).lower_bound = -flux
    
        # solve fba under these constraints
        solution = model.optimize()
    
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
    
    