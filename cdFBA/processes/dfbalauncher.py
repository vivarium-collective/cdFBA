import random
import pprint
import math

from process_bigraph import ProcessTypes, Process, Step, Composite
from process_bigraph.emitter import gather_emitter_results

from cdFBA.utils import model_from_file, get_single_dfba_spec, set_concentration, make_cdfba_composite, set_kinetics, get_objective_reaction
from cdFBA.processes.dfba import dFBA, UpdateEnvironment

from matplotlib import pyplot as plt

class EnvironmentMonitor(Step):
    """
    Monitors parameters in the shared environment and sends
    """
    config_schema = {}

    def __init__(self, config):
        super().__init__(self, config)

    def inputs(self):
        return {
            "thresholds": "map[threshold]",
            "shared_environment": "volumetric",
            "species": "any" #connect to "Species" store with all dFBAs
        }

    def outputs(self):
        return {
            "new_species":"map"
        }

    def update(self, inputs):

        to_add = {}
        to_remove = []

        for name, threshold in inputs["thresholds"].items():
            substrate = threshold["substrate"]
            if (threshold["upper"] is not None and inputs["shared_environment"][substrate] > threshold["upper"]) or (threshold["lower"] is not None and inputs["shared_environment"][substrate] < threshold["lower"]):
                if threshold["type"] == "add":
                    name = threshold["name"]
                    model = threshold["model"]
                    if threshold["parent"] is not None:
                        interval = inputs["species"][threshold["parent"]]["interval"]
                        config = inputs["species"][threshold["parent"]]["config"]
                    else:
                        raise ValueError("Please provide a parent species")
                    config["biomass_indentifier"] = get_objective_reaction(model_file=model)
                    config["model_file"] = model
                    config["name"] = name
                    spec = get_single_dfba_spec(model_file = model, name=name, config=config, interval=interval)
                    to_add[name] = spec
                if threshold["type"] == "remove":
                    model = threshold["model"]
                    to_remove.append(model)

        return {
            "manage_dfba": {
                "_add": to_add,
                "_remove": to_remove
            }
        }

def get_env_monitor_spec(interval):
    """Returns a specification dictionary for the environment monitor"""
    return {
        "_type": "process",
        "address": "local:EnvironmentMonitor",
        "config": {},
        "inputs": {
            "thresholds": ["thresholds"],
            "shared_environment": ["shared environment"],
            "species": ["Species"]
        },
        "outputs": {
            "manage_dfba": ["Species"]
        },
    }

def run_env_monitor(core):
    # BiGG model ids or the path name to the associated model file
    model_dict = {
        "E.coli": "iAF1260",
        # "S.flexneri": "iSFxv_1172"
    }
    # list exchange reactions
    exchanges = ["EX_glc__D_e", "EX_ac_e"]
    # set environment volume
    volume = 2
    # define a single dFBA model
    spec = make_cdfba_composite(model_dict, medium_type=None, exchanges=exchanges, volume=volume, interval=1.0)
    #create thresholds store
    spec["thresholds"] = {
        "Acetate": {
            "substrate": "Acetate",
            "range": {
                "upper": 20,
                "lower": None,
            },
            "model": "iAF1260",
            "parent": "E.coli",
            "name": "E.coli 2"
        }
    }
    spec["monitor"] = get_env_monitor_spec(interval=1.0)

    # Set reaction bounds
    spec["Species"]["E.coli"]["config"]["bounds"] = {
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
            "shared_environment": ["shared environment"],
            "global_time": ["global_time"]
        }
    }
    pprint.pprint(spec)

    # put it in a composite
    sim = Composite({
        "state": spec,
        # "emitter": {"mode": "all"}
    },
        core=core
    )
    pprint.pprint(sim.state)

    # run the simulation
    sim.run(40)
    results = gather_emitter_results(sim)[("emitter",)]

    # print the results
    timepoints = []
    for timepoint in results:
        time = timepoint.pop("global_time")
        timepoints.append(time)
        # dfba_spec = timepoint.pop(name1)
        print(f"TIME: {time}")
        print(f"STATE: {timepoint}")

    env = [timepoint["shared_environment"]["concentrations"] for timepoint in results]
    env_combined = {}
    for d in env:
        for key, value in d.items():
            if key not in env_combined:
                env_combined[key] = []
            env_combined[key].append(value)

    fig, ax = plt.subplots(dpi=300)
    for key, value in env_combined.items():
        # if not key == "glucose":
        #     continue
        ax.plot(timepoints, env_combined[key], label=key)
    plt.xlabel("Time")
    plt.ylabel("Substrate Concentration")
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    from cdFBA import register_types

    # create the core object
    core = ProcessTypes()
    # register data types
    core = register_types(core)
    # register all processes and steps
    core.register_process('dFBA', dFBA)
    core.register_process('UpdateEnvironment', UpdateEnvironment)
    core.register_process('EnvironmentMonitor', EnvironmentMonitor)

    run_env_monitor(core)