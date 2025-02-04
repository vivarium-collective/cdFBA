import random
import pprint
import math

from process_bigraph.composite import ProcessTypes
from process_bigraph import Process, Step, Composite

from cdFBA.utils import DFBAconfig, model_from_file, get_objective_reaction


from matplotlib import pyplot as plt


class DFBA(Process):
    """Performs single time-step of dynamic FBA

    Parameters:
    -----------
    model_file: string, math to cobra model file
    """
    config_schema = {
        "model_file": "any",
        "name": "any",
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
        "time_step": 0.1,
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

class Chemostat(Process):

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

class WaveFunction(Process):

    config_schema = {
        "substrate_params" : "map[map[float]]",
    }

    # substrate_params = {
    #     "substrate_name": {
    #         "amplitude": "float",
    #         "angular_frequency": "float",
    #         "base_concentration": "float",
    #         "phase_shift": "float"
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

class Injector(Process):

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

def run_environment(core):
    """This tests that the environment runs"""
    name1 = "E.coli"
    name2 = "S.flexneri"
    # define a single dFBA model
    spec = {
        "dfba": get_single_dfba_spec(model_file= "iAF1260", name=name1)
    }

    spec["dfba2"] = get_single_dfba_spec(model_file = "iSFxv_1172", name=name2)

    spec['shared environment'] = initial_environment(volume=2, species_list=[name1, name2])

    spec['dFBA Results'] = {name1:
        {
            "glucose": 0,
            "acetate": 0,
            spec['dfba']['config']['name']: 0,
        },
        name2:
        {
            "glucose": 0,
            "acetate": 0,
            spec['dfba2']['config']['name']: 0,
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
    # sim.update({
    #     "shared environment": {
    #         "glucose": 10,
    #         "acetate": 0,
    #         spec['dfba']['config']['biomass_identifier']: 0.1
    #         # "biomass": 0.1,
    #     }}
    #            )

    # get the results
    results = sim.gather_results()[('emitter',)]


    # print the results
    timepoints = []
    for timepoint in results:
        time = timepoint.pop('global_time')
        timepoints.append(time)
        dfba_spec = timepoint.pop('dfba')
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
        if not key == 'glucose':
            continue
        ax.plot(timepoints, env_combined[key], label=key)
    plt.xlabel('Time')
    plt.ylabel('Substrate Concentration')
    plt.legend()
    plt.tight_layout()
    plt.show()

def run_composite():
    pass


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
