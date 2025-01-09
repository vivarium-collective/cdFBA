import random
import pprint

from process_bigraph.composite import ProcessTypes
from process_bigraph import Process, Step, Composite

from cdFBA.utils import DFBAconfig, model_from_file


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
            "shared_environment": "map[float]", #initial conditions for time-step
            "current_update": "map[map[set_float]]",
        }

    def outputs(self):
        return {
             "dfba_update": "map[set_float]"
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
        state_update[self.config['biomass_identifier']] = biomass_growth_rate * current_biomass * interval

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
            "acetate": "EX_ac_e"}
    if bounds is None:
        bounds = {
            "EX_o2_e": {"lower": -2, "upper": None},
            "ATPM": {"lower": 1, "upper": 1}
        }

    if kinetics is None:
        kinetics = {
            "glucose": (0.5, 1),
            "acetate": (0.5, 2)}
    if biomass_identifier is None:
        biomass_identifier = DFBAconfig.get_objective_reaction(model=model)

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
        }
    }

class UpdateEnvironment(Step):
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
             "shared_environment": "map[float]",
             "species_updates": "map[map[set_float]]"
        }

    def outputs(self):
        return {
             "shared_environment": "map[float]"
        }

    def update(self, inputs):
        species_updates = inputs["species_updates"]
        shared_environment = inputs["shared_environment"]

        species_list = [species for species in species_updates]
        random.shuffle(species_list)

        update = shared_environment.copy()

        for species in species_list:
            for substrate_id in species_updates[species]:
                if (shared_environment[substrate_id] + species_updates[species][substrate_id]) > 0:
                    update[substrate_id] = species_updates[species][substrate_id]
                else:
                    update[substrate_id] = -shared_environment[substrate_id]

        # update = {}
        return {"shared_environment": update}

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
            "shared_environment": ["shared environment"]
        }
    }

def community_dfba_spec(
        species_list = [],
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

def run_dfba_alone(core):

    model_file = "textbook"
    config = dfba_config(model_file=model_file)
    dfba = DFBA(config, core=core)

    initial_state = dfba.initial_state()
    inputs = {"shared_environment": initial_state}

    results = dfba.update(inputs, 1)

    print(results)

def run_dfba(core):
    name = "E.coli"
    # define a single dFBA model
    spec = {
        "dfba": get_single_dfba_spec(name=name)
    }

    # TODO -- more automatic way to get initial environment
    spec['shared environment'] = {
        "glucose": 10,
        "acetate": 0,
        spec['dfba']['config']['biomass_identifier']: 0.1
        # "biomass": 0.1,
    }

    spec['dFBA Results'] = {
        name:
        {
            "glucose": 0,
            "acetate": 0,
            spec['dfba']['config']['biomass_identifier']: 0
        }
    }

    # put it in a composite
    sim = Composite({
        "state": spec,
        "emitter": {'mode': 'all'}},
        core=core
    )
    print(spec)
    # run the simulation
    sim.run(10)

    # get the results
    results = sim.gather_results()[('emitter',)]

    # print the results
    for timepoint in results:
        time = timepoint.pop('global_time')
        dfba_spec = timepoint.pop('dfba')
        print(f'TIME: {time}')
        print(f'STATE: {timepoint}')


def run_environment(core):
    """This tests that the environment runs"""
    name = "E.coli"
    # define a single dFBA model
    spec = {
        "dfba": get_single_dfba_spec(model_file= "iAF1260", name=name)
    }

    spec['shared environment'] = {
        "glucose": 100,
        "acetate": 5,
        spec['dfba']['config']['biomass_identifier']: 0.5
        # "biomass": 0.1,
    }

    spec['dFBA Results'] = {name:
        {
        "glucose": 0,
        "acetate": 0,
        spec['dfba']['config']['biomass_identifier']: 0
        }
    }

    spec['update environment'] = environment_spec()
    # put it in a composite
    sim = Composite({
        "state": spec,
        "emitter": {'mode': 'all'}},
        core=core
    )
    pprint.pprint(spec)

    # run the simulation
    sim.run(10)
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
    for timepoint in results:
        time = timepoint.pop('global_time')
        dfba_spec = timepoint.pop('dfba')
        print(f'TIME: {time}')
        print(f'STATE: {timepoint}')
    pass


def run_composite():
    pass


if __name__ == "__main__":
    from cdFBA import register_types

    # create a core
    core = ProcessTypes()
    core = register_types(core)

    core.register_process('DFBA', DFBA)
    core.register_process('UpdateEnvironment', UpdateEnvironment)

    # print(get_single_dfba_spec())
    # test_dfba_alone(core)
    # test_dfba(core)
    run_environment(core)
    # test_composite()
