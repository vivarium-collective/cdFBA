import random
import pprint
import math

from process_bigraph.composite import ProcessTypes
from process_bigraph import Process, Step, Composite

from cdFBA.utils import model_from_file, get_objective_reaction, get_injector_spec, get_wave_spec, get_chemo_spec
from cdFBA.utils import get_single_dfba_spec, dfba_config, environment_spec, initial_environment

from matplotlib import pyplot as plt


class DFBA(Process):
    """Performs single time-step of dynamic FBA

    Parameters:
    -----------
    model_file: string, math to cobra model file
    """
    config_schema = {
        "model_file": {
            '_type': 'string',
            '_default': 'iAF1260',
        },
        "name": 'string',
        "kinetics": "any",
        "reaction_map": "any",
        "biomass_identifier": "any",
        "bounds": "maybe[map[bounds]]",
    }

    def __init__(self, config, core):
        super().__init__(config, core)

        # TODO -- load model here rather than passing it in
        self.model = model_from_file(self.config['model_file'])


        if self.config["bounds"] is not None:
            if len(self.config["bounds"]) != 0:
                for reaction_id, bounds in self.config["bounds"].items():
                    if bounds["lower"] is not None:
                        self.model.reactions.get_by_id(reaction_id).lower_bound = bounds["lower"]
                    if bounds["upper"] is not None:
                        self.model.reactions.get_by_id(reaction_id).upper_bound = bounds["upper"]

    # def initial_state(self):
    #     # TODO -- get the initial state from the load model, self.model
    #     return {
    #         "dfba_update": {
    #             'Biomass_Ecoli_core': 0,
    #             'acetate': 0,
    #             'glucose': 0}}

    def inputs(self):
        return {
            "shared_environment": "volumetric", #initial conditions for time-step
            "current_update": "map[map[set_float]]",
        }

    def outputs(self):
        return {
             "dfba_update": "map[set_float]"
        }


    def update(self, inputs, interval):
        current_state = {key:inputs["shared_environment"]["counts"][key] for key, value in self.config["reaction_map"].items()}
        current_state[self.config["name"]] = inputs["shared_environment"]["counts"][self.config["name"]]
        state_update = current_state.copy()

        for substrate_id, reaction_id in self.config["reaction_map"].items():
            Km, Vmax = self.config["kinetics"][substrate_id]
            substrate_concentration = inputs["shared_environment"]["concentrations"][substrate_id]

            # calculate michaelis-menten flux
            flux = Vmax * substrate_concentration / (Km + substrate_concentration)

            # use the flux to constrain fba
            self.model.reactions.get_by_id(reaction_id).lower_bound = -flux

        # solve fba under these constraints
        solution = self.model.optimize()

        # gather the results
        ## update biomass
        biomass_growth_rate = solution.fluxes[self.config["biomass_identifier"]]
        current_biomass = current_state[self.config["name"]]
        state_update[self.config['name']] = biomass_growth_rate * current_biomass * interval

        ## update substrates
        for substrate_id, reaction_id in self.config["reaction_map"].items():
            flux = solution.fluxes[reaction_id]
            state_update[substrate_id] = (flux * current_biomass * interval)

        return {"dfba_update": state_update}

class UpdateEnvironment(Step): #TODO =:
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
             "shared_environment": "volumetric",
             "species_updates": "map[map[set_float]]" #TODO add a species_update type
        }

    def outputs(self):
        return {
            "shared_environment": "map[float]",
        }

    def update(self, inputs):
        species_updates = inputs["species_updates"]
        shared_environment = inputs["shared_environment"]["counts"]
        env_volume = inputs["shared_environment"]["volume"]

        species_list = [species for species in species_updates]
        random.shuffle(species_list)

        update = shared_environment.copy()

        for species in species_list:
            for substrate_id in species_updates[species]:
                if (shared_environment[substrate_id] + species_updates[species][substrate_id]) > 0:
                    update[substrate_id] = species_updates[species][substrate_id]
                else:
                    update[substrate_id] = -shared_environment[substrate_id]

        return {
            "shared_environment": {'counts': update}
        }

class Chemostat(Process):
    """The Chemostat process maintains the concentration of given substrates at a fixed value at each time-step
    """
    config_schema = {
        "substrate_concentrations" : "map[float]",
    }

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "shared_environment": "volumetric",
            "global_time": "float"
        }

    def outputs(self):
        return {
            "shared_environment": "map[float]"
        }

    def update(self, inputs, interval):
        shared_environment = inputs["shared_environment"]["counts"]

        update = {}

        for substrate, values in self.config["substrate_concentrations"].items():
            update[substrate] = (values * inputs["shared_environment"]["volume"]) - shared_environment[substrate]

        return {
            "shared_environment": {'counts': update}
        }

class WaveFunction(Process):
    """The WaveFunction process maintains the concentration of given substrates based on a wave-function
    """
    config_schema = {
        "substrate_params" : "map[map[float]]",
    }

    """
    substrate_params = {
        "substrate_name": {
            "amplitude": "float",
            "angular_frequency": "float",
            "base_concentration": "float",
            "phase_shift": "float"
        }
    }
    
    amplitude : float, amplitude of wave function
    angular_frequency : float, angular frequency of wave function
    base_concentration : float, base concentration substrate
    phase_shift : float, phase shift of wave function (when wave starts)
    """

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "shared_environment": "volumetric",
            "global_time": "float"
        }

    def outputs(self):
        return {
            "shared_environment": "map[float]"
        }

    def update(self, inputs, interval):
        shared_environment = inputs["shared_environment"]["counts"]
        t = inputs["global_time"]
        update = {}
        for substrate in self.config["substrate_params"]:
            A = self.config["substrate_params"][substrate]["amplitude"]
            w = self.config["substrate_params"][substrate]["angular_frequency"]
            B = self.config["substrate_params"][substrate]["base_concentration"]
            phi = self.config["substrate_params"][substrate]["phase_shift"]

            current_count = (A*math.sin(w*t+phi) + B) * inputs["shared_environment"]["volume"]

            update[substrate] = (current_count - shared_environment[substrate])

        return {
            "shared_environment": {'counts': update}
        }

class Injector(Process):
    """The Injector process injects a given amount of a given substrate at regular intervals into the shared environment
    """
    config_schema = {
        "injection_params" : "map[map[float]]",
    }

    # injection_params = {
    #     "substrate_name" : {
    #         "amount": "float",
    #         "interval": "float",
    #     }
    # }

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "shared_environment": "volumetric",
            "global_time": "float"
        }

    def outputs(self):
        return {
            "shared_environment": "map[float]"
        }

    def update(self, inputs, interval):
        shared_environment = inputs["shared_environment"]["counts"]
        t = inputs["global_time"]
        update = {}
        for substrate in self.config["injection_params"]:
            if ((t % self.config["injection_params"][substrate]["interval"]) == 0) & (t!=0.0):
                update[substrate] = self.config["injection_params"][substrate]["amount"]
        return {
            "shared_environment": {'counts': update}
        }

def run_environment(core):
    """This tests that the environment runs"""
    name1 = "E.coli"
    name2 = "S.flexneri"
    # define a single dFBA model
    spec = {
        name1: get_single_dfba_spec(model_file= "iAF1260", name=name1)
    }

    spec[name2] = get_single_dfba_spec(model_file = "iSFxv_1172", name=name2)

    spec['shared environment'] = initial_environment(volume=2, species_list=[name1, name2])

    spec['dFBA Results'] = {name1:
        {
            "glucose": 0,
            "acetate": 0,
            spec[name1]['config']['name']: 0,
        },
        name2:
        {
            "glucose": 0,
            "acetate": 0,
            spec[name2]['config']['name']: 0,
        }
    }

    spec['update environment'] = environment_spec()

    injector_config = {
        "injection_params": {
            "glucose": {
                "amount": 80,
                "interval": 5,
            }
        }
    }

    spec['environment dynamics'] = get_injector_spec(config=injector_config)

    pprint.pprint(spec)

    # put it in a composite
    sim = Composite({
        "state": spec,
        "emitter": {'mode': 'all'}},
        core=core
    )

    # run the simulation
    sim.run(40)
    results = sim.gather_results()[('emitter',)]

    # print the results
    timepoints = []
    for timepoint in results:
        time = timepoint.pop('global_time')
        timepoints.append(time)
        dfba_spec = timepoint.pop(name1)
        print(f'TIME: {time}')
        print(f'STATE: {timepoint}')

    env = [timepoint['shared environment']['concentrations'] for timepoint in results]
    env_combined = {}
    for d in env:
        for key, value in d.items():
            if key not in env_combined:
                env_combined[key] = []
            env_combined[key].append(value)


    fig, ax = plt.subplots(dpi=300)
    for key, value in env_combined.items():
        # if not key == 'glucose':
        #     continue
        ax.plot(timepoints, env_combined[key], label=key)
    plt.xlabel('Time')
    plt.ylabel('Substrate Concentration')
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    from cdFBA import register_types

    # create a core
    core = ProcessTypes()
    core = register_types(core)

    core.register_process('DFBA', DFBA)
    core.register_process('UpdateEnvironment', UpdateEnvironment)
    core.register_process('Chemostat', Chemostat)
    core.register_process('WaveFunction', WaveFunction)
    core.register_process('Injector', Injector)

    # print(get_single_dfba_spec())
    # test_dfba_alone(core)
    # test_dfba(core)
    run_environment(core)
    # test_composite()
