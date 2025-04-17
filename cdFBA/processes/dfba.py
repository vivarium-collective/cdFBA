import random
import pprint
import math
import pytest

from process_bigraph import Process, Step, Composite, ProcessTypes
from process_bigraph.emitter import gather_emitter_results

from cdFBA.utils import SHARED_ENVIRONMENT
from cdFBA.utils import model_from_file, get_injector_spec, get_wave_spec, get_static_spec, set_concentration
from cdFBA.utils import  make_cdfba_composite, set_kinetics, get_objective_reaction

from matplotlib import pyplot as plt


class dFBA(Process):
    """Performs single time-step of dynamic FBA

    Config Parameters:
    -----------
    model_file: str, big model ID or path to cobra model file
    name: string, name of process (usually species/strain name)
    kinetics: dict, dictionary of tuples with kinetic parameters (km, Vmax)
    reaction_map: dict, maps substrate names to reaction IDs
    bounds: dict, maps reaction IDs to a bounds dictionary
    """
    config_schema = {
        "model_file": {
            "_type": "any",
            "_default": "iAF1260",
        },
        "name": "string",
        "kinetics": "any",
        "reaction_map": "any",
        "bounds": "maybe[map[bounds]]",
        "changes": "dfba_changes"
    }
    #TODO -- Add "changes" to the config that allows us to change the model
    def __init__(self, config, core):
        super().__init__(config, core)

        self.model = model_from_file(self.config["model_file"])
        self.biomass_identifier = get_objective_reaction(self.model)

        if self.config["bounds"] is not None:
            if len(self.config["bounds"]) != 0:
                for reaction_id, bounds in self.config["bounds"].items():
                    if bounds["lower"] is not None:
                        self.model.reactions.get_by_id(reaction_id).lower_bound = bounds["lower"]
                    if bounds["upper"] is not None:
                        self.model.reactions.get_by_id(reaction_id).upper_bound = bounds["upper"]

        if self.config["changes"] is not None:
            if len(self.config["changes"]["gene_knockout"]) > 0:
                for gene in self.config["changes"]["gene_knockout"]:
                    self.model.genes.get_by_id(gene).knock_out()
            if len(self.config["changes"]["reaction_knockout"]) > 0:
                for reaction in self.config["changes"]["reaction_knockout"]:
                    self.model.reactions.get_by_id(reaction).knock_out()
            if len(self.config["changes"]["bounds"]) > 0:
                for reaction_id, bounds in self.config["changes"]["bounds"].items():
                    if bounds["lower"] is not None:
                        self.model.reactions.get_by_id(reaction_id).lower_bound = bounds["lower"]
                    if bounds["upper"] is not None:
                        self.model.reactions.get_by_id(reaction_id).upper_bound = bounds["upper"]
            if len(self.config["changes"]["kinetics"]) > 0:
                self.config["kinetics"].update(self.config["changes"]["kinetics"])

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

            # calculate Michaelis-Menten flux
            flux = Vmax * substrate_concentration / (Km + substrate_concentration)

            # use the flux to constrain fba
            self.model.reactions.get_by_id(reaction_id).lower_bound = -flux

        # solve fba under these constraints
        solution = self.model.optimize()

        # gather the results
        ## update biomass
        biomass_growth_rate = solution.fluxes[self.biomass_identifier]
        current_biomass = current_state[self.config["name"]]
        state_update[self.config["name"]] = biomass_growth_rate * current_biomass * interval

        ## update substrates
        for substrate_id, reaction_id in self.config["reaction_map"].items():
            flux = solution.fluxes[reaction_id]
            state_update[substrate_id] = (flux * current_biomass * interval)

        return {"dfba_update": state_update}

class UpdateEnvironment(Step):
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "shared_environment": "volumetric",
            "species_updates": "map[map[set_float]]",
        }

    def outputs(self):
        return {
            "counts": "map[float]",
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
            "counts": update
        }

class StaticConcentration(Process):
    """The StaticConcentration process maintains the concentration of given substrates at a fixed value at each time-step
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
            "shared_environment": {"counts": update}
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
            if current_count > 0:
                update[substrate] = (current_count - shared_environment[substrate])
            else:
                update[substrate] = 0

        return {
            "shared_environment": {"counts": update}
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
            if (((t+1) % self.config["injection_params"][substrate]["interval"]) == 0) & (t!=0.0):
                update[substrate] = self.config["injection_params"][substrate]["amount"]
        return {
            "shared_environment": {"counts": update}
        }

#=======
# TESTS
#=======

def get_test_spec():
    # BiGG model ids or the path name to the associated model file
    model_dict = {
        "E.coli": "iAF1260",
        "S.flexneri": "iSFxv_1172"
    }
    # list exchange reactions
    exchanges = ["EX_glc__D_e", "EX_ac_e"]
    # set environment volume
    volume = 2
    # define a single dFBA model
    spec = make_cdfba_composite(model_dict, medium_type=None, exchanges=exchanges, volume=volume, interval=1.0)

    # Set reaction bounds
    spec["Species"]["E.coli"]["config"]["bounds"] = {
        "EX_o2_e": {"lower": -2, "upper": None},
        "ATPM": {"lower": 1, "upper": 1}
    }
    spec["Species"]["S.flexneri"]["config"]["bounds"] = {
        "EX_o2_e": {"lower": -2, "upper": None},
        "ATPM": {"lower": 1, "upper": 1}
    }

    # set external substrate concentrations
    concentrations = {
        "Acetate": 0,
        "D-Glucose": 40
    }
    set_concentration(spec, concentrations)

    # set kinetics
    kinetics = {
        "D-Glucose": (0.02, 15),
        "Acetate": (0.5, 7)
    }
    for species in model_dict.keys():
        set_kinetics(species, spec, kinetics)

    # set emitter specs
    spec["emitter"] = {
        "_type": "step",
        "address": "local:ram-emitter",
        "config": {
            "emit": {
                "shared_environment": "any",
                "global_time": "any",
            }
        },
        "inputs": {
            "shared_environment": [SHARED_ENVIRONMENT],
            "global_time": ["global_time"]
        }
    }
    return spec

@pytest.fixture
def core():
    from cdFBA import register_types
    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    # register all processes and steps
    core.register_process("dFBA", dFBA)
    core.register_process("UpdateEnvironment", UpdateEnvironment)
    core.register_process("StaticConcentration", StaticConcentration)
    core.register_process("WaveFunction", WaveFunction)
    core.register_process("Injector", Injector)

    return core

def test_environment(core):
    """This tests that the environment runs"""
    spec = get_test_spec()
    # put it in a composite
    sim = Composite({
        "state": spec,
        # "emitter": {"mode": "all"}
    },
        core=core
    )

    # run the simulation
    sim.run(20)
    results = gather_emitter_results(sim)[("emitter",)]

    assert len(results) == 21
    assert results[2]["shared_environment"]["counts"]["D-Glucose"] > results[4]["shared_environment"]["counts"]["D-Glucose"]
    assert results[4]["shared_environment"]["counts"]["E.coli"] > results[2]["shared_environment"]["counts"]["E.coli"]

def test_static_concentration(core):
    spec = get_test_spec()
    StaticConcentration_config = {
        "substrate_concentrations": {
            "D-Glucose": 40
        }
    }
    spec['StaticConcentration'] = get_static_spec(StaticConcentration_config, interval=1.0)
    # put it in a composite
    sim = Composite({
        "state": spec,
        # "emitter": {"mode": "all"}
    },
        core=core
    )

    # run the simulation
    sim.run(10)
    results = gather_emitter_results(sim)[("emitter",)]

    assert len(results) == 11
    assert results[2]["shared_environment"]["counts"]["D-Glucose"] == results[4]["shared_environment"]["counts"][
        "D-Glucose"]
    assert results[4]["shared_environment"]["counts"]["E.coli"] > results[2]["shared_environment"]["counts"]["E.coli"]

def test_injector(core):
    spec = get_test_spec()
    injector_config = {
        "injection_params": {
            "D-Glucose": {
                "amount": 80,
                "interval": 5,
            }
        }
    }
    spec['Injector'] = get_injector_spec(injector_config, interval=1.0)
    # put it in a composite
    sim = Composite({
        "state": spec,
        # "emitter": {"mode": "all"}
    },
        core=core
    )

    # run the simulation
    sim.run(20)
    results = gather_emitter_results(sim)[("emitter",)]

    assert len(results) == 21
    assert results[4]["shared_environment"]["counts"]["E.coli"] > results[2]["shared_environment"]["counts"]["E.coli"]
    assert results[10]["shared_environment"]["counts"]["D-Glucose"]== results[20]["shared_environment"]["counts"]["D-Glucose"]

if __name__ == "__main__":
    from cdFBA import register_types

    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    # register all processes and steps
    core.register_process("dFBA", dFBA)
    core.register_process("UpdateEnvironment", UpdateEnvironment)
    core.register_process("StaticConcentration", StaticConcentration)
    core.register_process("WaveFunction", WaveFunction)
    core.register_process("Injector", Injector)

    test_environment(core)
    test_static_concentration(core)
    test_injector(core)
