import random
from math import isclose, sin
from itertools import cycle

from process_bigraph import Process, Step

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
    """The single writer of the shared environment's counts.

    Each step it (1) applies the dFBA consumption/production deltas, clamped so no
    pool goes negative, then (2) enforces any held-concentration targets as boundary
    conditions: a substrate named in ``concentration_targets`` is forced to
    ``target * volume`` after the dynamics, which replenishes whatever was consumed
    (a chemostat). Making this the ONLY process that writes ``shared_environment``
    is what makes a hold genuine: the setter processes (StaticConcentration /
    WaveFunction) merely DECLARE targets into a separate store, so there is no
    write-write cycle between two steps on ``shared_environment`` (which is what let
    the consumption land last and the hold drift away, hidden by emitter timing).
    """
    config_schema = {}

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "shared_environment": "volumetric",
            "species_updates": "map[map[overwrite[float]]]",
            "concentration_targets": "map[overwrite[float]]",
        }

    def outputs(self):
        return {
            "shared_environment": "volumetric",
        }

    def update(self, inputs):
        species_updates = inputs["species_updates"]
        shared_environment = inputs["shared_environment"]["counts"]
        env_volume = inputs["shared_environment"]["volume"]
        targets = inputs.get("concentration_targets") or {}

        species_list = [species for species in species_updates]
        random.shuffle(species_list)

        update = {substrate: 0 for substrate in shared_environment}

        # 1) apply dFBA consumption/production. Species draw in random order,
        # first-come-first-served: each is clamped against what is actually left
        # after the preceding species, so a substrate depletes to exactly zero and
        # never goes negative.
        for species in species_list:
            for substrate_id, delta in species_updates[species].items():
                remaining = shared_environment[substrate_id] + update[substrate_id]
                if remaining + delta < 0:
                    update[substrate_id] = -shared_environment[substrate_id]
                else:
                    update[substrate_id] += delta

        # 2) enforce held concentrations as boundary conditions AFTER the dynamics:
        # force each held substrate's settled count to target * volume, overriding
        # this step's consumption (the chemostat inflow that replenishes it). A
        # non-positive target means "not held this step" (e.g. a wave dipping <= 0),
        # so the substrate is left to its consumption dynamics.
        for substrate_id, target_conc in targets.items():
            if substrate_id in shared_environment and target_conc > 0:
                target_count = target_conc * env_volume
                update[substrate_id] = target_count - shared_environment[substrate_id]

        return {
            "shared_environment": {
                "counts": update
            }
        }

class StaticConcentration(Process):
    """Declares fixed concentration targets for given substrates.

    A pure target-declarer: it writes each held substrate's target concentration to
    the shared ``concentration_targets`` store and never touches ``shared_environment``.
    UpdateEnvironment (the sole writer of the environment) enforces those targets after
    the tick's consumption, which is what makes the hold genuine and keeps this process
    out of any write-write cycle on the environment counts.
    """
    config_schema = {
        "substrate_concentrations" : "map[float]",
    }

    def __init__(self, config, core):
        super().__init__(config, core)

    def inputs(self):
        return {
            "global_time": "float",
        }

    def outputs(self):
        return {
            "concentration_targets": "map[overwrite[float]]",
        }

    def update(self, inputs, interval):
        return {
            "concentration_targets": dict(self.config["substrate_concentrations"]),
        }

class WaveFunction(Process):
    """Declares a time-varying (wave) concentration target for given substrates.

    Like StaticConcentration, a pure target-declarer into ``concentration_targets``
    (never touches ``shared_environment``). The target concentration at time t is
    ``A*sin(w*t + phi) + B``; a non-positive target is written as 0, which
    UpdateEnvironment reads as "not held this step" so the substrate follows its
    consumption dynamics rather than being pinned to a negative concentration.
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
            "global_time": "float"
        }

    def outputs(self):
        return {
            "concentration_targets": "map[overwrite[float]]",
        }

    def update(self, inputs, interval):
        t = inputs["global_time"]
        targets = {}
        for substrate in self.config["substrate_params"]:
            A = self.config["substrate_params"][substrate]["amplitude"]
            w = self.config["substrate_params"][substrate]["angular_frequency"]
            B = self.config["substrate_params"][substrate]["base_concentration"]
            phi = self.config["substrate_params"][substrate]["phase_shift"]

            target_concentration = A * sin(w * t + phi) + B
            targets[substrate] = target_concentration if target_concentration > 0 else 0.0

        return {
            "concentration_targets": targets,
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
