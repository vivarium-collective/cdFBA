import random
from math import isclose, sin
from itertools import cycle

from process_bigraph import Process, Step

from viva_cobra import mm_uptake_rates
from viva_cobra.processes.fba import fba_deltas

from cdFBA.utils import model_from_file, get_objective_reaction


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
            "_type": "string",
            "_default": "iAF1260",
        },
        "name": "string",
        "kinetics": "map",
        "reaction_map": "map",
        "bounds": "maybe[map[bounds]]",
        "changes": "dfba_changes",
        "medium": "maybe[map]"
    }
    #TODO -- add ability to change objective reaction
    def __init__(self, config, core):
        super().__init__(config, core)

        self.model = model_from_file(self.config["model_file"])
        self.biomass_identifier = get_objective_reaction(self.model)

        # `maybe` config fields that are unset are absent from config entirely,
        # so read every optional field defensively.
        medium = self.config.get("medium") or {}
        if len(medium) > 0:
            self.model.medium = medium

        self._apply_bounds(self.config.get("bounds") or {})

        changes = self.config.get("changes") or {}
        for gene in changes.get("gene_knockout") or []:
            self.model.genes.get_by_id(gene).knock_out()
        for reaction in changes.get("reaction_knockout") or []:
            self.model.reactions.get_by_id(reaction).knock_out()
        self._apply_bounds(changes.get("bounds") or {})
        self.config["kinetics"].update(changes.get("kinetics") or {})

    def _apply_bounds(self, bounds_map):
        """Apply a {reaction_id: {"lower": .., "upper": ..}} map to the model."""
        for reaction_id, bounds in bounds_map.items():
            if bounds.get("lower") is not None:
                self.model.reactions.get_by_id(reaction_id).lower_bound = bounds["lower"]
            if bounds.get("upper") is not None:
                self.model.reactions.get_by_id(reaction_id).upper_bound = bounds["upper"]

    def inputs(self):
        return {
            "shared_environment": "volumetric", #initial conditions for time-step
            "current_update": "map[map[overwrite[float]]]",
        }

    def outputs(self):
        return {
             "dfba_update": "map[overwrite[float]]"
        }

    def update(self, inputs, interval):
        # Reuse viva-cobra's FBA process instead of solving COBRApy here: Michaelis-Menten
        # uptake, then viva-cobra's `fba_deltas` (the same LP solve, additionally supply-capped
        # so no pool goes negative, and starvation-safe). The community coupling — the shared
        # environment and the per-species `dfba_update` — stays this workspace's glue.
        name = self.config["name"]
        reaction_map = self.config["reaction_map"]          # {substrate: exchange reaction id}
        resources = list(reaction_map)
        env = inputs["shared_environment"]
        substrates = {s: env["concentrations"][s] for s in resources}
        biomass = env["counts"][name]

        uptakes = mm_uptake_rates(substrates, self.config["kinetics"], resources)
        deltas = fba_deltas(self.model, reaction_map, uptakes, substrates, biomass, interval)

        # viva-cobra returns {"substrates": {resource: Δ}, "biomass": Δ}; cdFBA's dfba_update
        # carries the per-species biomass delta under the species name alongside the substrates.
        state_update = dict(deltas["substrates"])
        state_update[name] = deltas["biomass"]
        return {"dfba_update": state_update}

class UpdateEnvironment(Step):
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "shared_environment": "volumetric",
            "species_updates": "map[map[overwrite[float]]]",
        }

    def outputs(self):
        return {
            "shared_environment": "volumetric",
        }

    def update(self, inputs):
        species_updates = inputs["species_updates"]
        shared_environment = inputs["shared_environment"]["counts"]
        env_volume = inputs["shared_environment"]["volume"]

        species_list = [species for species in species_updates]
        random.shuffle(species_list)

        update = {substrate: 0 for substrate in shared_environment}

        # species draw in random order, first-come-first-served: each one is
        # clamped against what is actually left after the preceding species,
        # so a substrate is depleted to exactly zero and never goes negative.
        for species in species_list:
            for substrate_id, delta in species_updates[species].items():
                remaining = shared_environment[substrate_id] + update[substrate_id]
                if remaining + delta < 0:
                    update[substrate_id] = -shared_environment[substrate_id]
                else:
                    update[substrate_id] += delta

        return {
            "shared_environment": {
                "counts": update
            }
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
            "shared_environment": "volumetric",
        }

    def update(self, inputs, interval):
        shared_environment = inputs["shared_environment"]["counts"]

        update = {}

        for substrate, values in self.config["substrate_concentrations"].items():
            update[substrate] = (values * inputs["shared_environment"]["volume"]) - shared_environment[substrate]

        return {
            "shared_environment": {
                "counts": update
            }
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
            "shared_environment": "volumetric",
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

            current_count = (A*sin(w*t+phi) + B) * inputs["shared_environment"]["volume"]
            if current_count > 0:
                update[substrate] = (current_count - shared_environment[substrate])
            else:
                update[substrate] = 0

        return {
            "shared_environment": {
                "counts": update
            }
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
    #         "interval": "float"/"list[float]",
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
            "shared_environment": "volumetric",
        }

    def update(self, inputs, interval):
        tol = 1e-6
        shared_environment = inputs["shared_environment"]["counts"]
        t = inputs["global_time"]
        update = {}

        for substrate in self.config["injection_params"]:
            intervals = self.config["injection_params"][substrate]["interval"]
            if isinstance(intervals, (float,int)):
                intervals = [intervals]

            total = 0.0
            for step in cycle(intervals):
                total += step
                if total > t+tol:
                    break
                if isclose(t, total, abs_tol=tol) and t != 0.0:
                    update[substrate] = self.config["injection_params"][substrate]["amount"]
        return {
            "shared_environment": {
                "counts": update
            }
        }
