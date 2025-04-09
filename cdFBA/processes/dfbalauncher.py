import random
import pprint
import math

from process_bigraph import ProcessTypes, Process, Step, Composite
from process_bigraph.emitter import gather_emitter_results

from cdFBA.utils import SHARED_ENVIRONMENT, SPECIES_STORE, THRESHOLDS, DFBA_RESULTS
from cdFBA.utils import model_from_file, get_single_dfba_spec, set_concentration, make_cdfba_composite, set_kinetics, get_objective_reaction
from cdFBA.processes.dfba import dFBA, UpdateEnvironment

from matplotlib import pyplot as plt

class EnvironmentMonitor(Step):
    """
    Monitors parameters in the shared environment and sends
    """
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "thresholds": "map[threshold]",
            "shared_environment": "volumetric",
            "species": "any", #connect to SPECIES_STORE store with all dFBAs
            "dfba_results": "any",
        }

    def outputs(self):
        return {
            "new_species": "map",
            "concentrations": "map",
            "counts": "map",
            "dfba_results": "any",
            "mass_removal": "map[float]"
        }

    def update(self, inputs):

        to_add = {}
        to_remove = []

        add_concentrations = {}
        add_counts = {}

        remove_concentrations = []
        remove_counts = []

        add_dfba_updates = {}
        remove_dfba_updates = []

        mass_updates = {}

        for threshold in inputs["thresholds"].values():
            substrate = threshold["substrate"]
            if ((isinstance(threshold["range"]["upper"], (float, int))
                and inputs["shared_environment"]["concentrations"][substrate] > threshold["range"]["upper"])
                    or (isinstance(threshold["range"]["lower"], (float, int))
                        and inputs["shared_environment"]["concentrations"][substrate] < threshold["range"]["lower"])):
                name = threshold["name"]
                parent = threshold["parent"]
                mass = threshold["mass"]
                if threshold["type"] == "add":
                    if not name in inputs["species"].keys():
                        interval = inputs["species"][parent]["interval"]
                        config = inputs["species"][parent]["config"]
                        config["name"] = name
                        config["changes"] = threshold["changes"]
                        spec = get_single_dfba_spec(model_file=config["model_file"], name=threshold["name"], config=config, interval=interval)
                        to_add[name] = spec
                        add_counts[name] = inputs["thresholds"][threshold]["counts"]
                        add_concentrations[name] = inputs["thresholds"][threshold]["mass"]/inputs["shared_environment"]["volume"]
                        environment_substrates = [substrate for substrate in inputs["dfba_results"][parent].keys() if substrate != parent]
                        environment_substrates.append(name)
                        add_dfba_updates[name] = {substrate: 0 for substrate in environment_substrates}
                        add_dfba_updates[name][substrate] = mass
                        mass_updates[parent] = {parent: -mass}

                if threshold["type"] == "remove":
                    to_remove.append(name)
                    remove_counts.append(name)
                    remove_concentrations.append(name)
                    remove_dfba_updates.append(name)

        if to_add:
            import ipdb; ipdb.set_trace()
        return {
            "new_species": {
                '_add': to_add,
                '_remove': to_remove
            },
            "concentrations": {
                '_add': add_concentrations,
                '_remove': remove_concentrations,
            },
            "counts": {
                '_add': add_counts,
                '_remove': remove_counts,
            },
            "dfba_results": {
                '_add': add_dfba_updates,
                '_remove': remove_dfba_updates,
            },
            "mass_removal": {"counts": mass_updates}
        }

def get_env_monitor_spec(interval):
    """Returns a specification dictionary for the environment monitor"""
    return {
        "_type": "process",
        "address": "local:EnvironmentMonitor",
        "config": {},
        "inputs": {
            "thresholds": [THRESHOLDS],
            "shared_environment": [SHARED_ENVIRONMENT],
            "species": [SPECIES_STORE],
            "dfba_results": [DFBA_RESULTS],
        },
        "outputs": {
            "new_species": [SPECIES_STORE],
            "concentrations": [SHARED_ENVIRONMENT, "concentrations"],
            "counts": [SHARED_ENVIRONMENT, "counts"],
            "dfba_results": [DFBA_RESULTS],
            "mass_removal": [SHARED_ENVIRONMENT],
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
    spec[THRESHOLDS] = {
        "Acetate": {
            "type": "add",
            "substrate": "Acetate",
            "range": {
                "upper": 20,
                "lower": None,
            },
            "parent": "E.coli",
            "name": "E.coli 2",
            "changes":{
                "gene_knockout": [],
                "reaction_knockout": [],
            }
        }
    }
    spec["monitor"] = get_env_monitor_spec(interval=1.0)

    # Set reaction bounds
    spec[SPECIES_STORE]["E.coli"]["config"]["bounds"] = {
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

    pprint.pprint(sim.state)

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
