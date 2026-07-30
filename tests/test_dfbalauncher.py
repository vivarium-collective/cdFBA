"""Tests for the EnvironmentMonitor step, which spawns and removes species."""
import pytest

from cdFBA.processes.dfbalauncher import EnvironmentMonitor, get_env_monitor_spec
from cdFBA.utils import SHARED_ENVIRONMENT, SPECIES_STORE, THRESHOLDS, DFBA_RESULTS
from tests.conftest import build_spec, env, run_spec, series

NO_CHANGES = {
    "gene_knockout": [],
    "reaction_knockout": [],
    "bounds": {},
    "kinetics": {},
}


def threshold(upper=None, lower=None, kind="add", substrate="Acetate",
              parent="E.coli", name="E.coli 2", mass=0.25):
    """Build a threshold value, omitting unset bounds as the schema does."""
    bounds = {}
    if upper is not None:
        bounds["upper"] = upper
    if lower is not None:
        bounds["lower"] = lower
    return {
        "type": kind,
        "substrate": substrate,
        "range": bounds,
        "parent": parent,
        "name": name,
        "changes": dict(NO_CHANGES),
        "mass": mass,
    }


def monitor_inputs(thresholds, concentrations, species=("E.coli",)):
    counts = dict(concentrations)
    return {
        "thresholds": thresholds,
        "shared_environment": {
            "volume": 1.0, "counts": counts, "concentrations": dict(concentrations)},
        "species": {
            name: {
                "interval": 1.0,
                "config": {"model_file": "textbook", "name": name, "changes": dict(NO_CHANGES)},
            } for name in species
        },
        "dfba_results": {
            name: {"D-Glucose": 0.0, "Acetate": 0.0, name: 0.0} for name in species
        },
    }


@pytest.fixture
def monitor(core):
    return EnvironmentMonitor({}, core)


def test_ports(monitor):
    assert set(monitor.inputs()) == {
        "thresholds", "shared_environment", "species", "dfba_results"}
    assert set(monitor.outputs()) == {"new_species", "counts", "dfba_results"}


def test_no_change_when_threshold_not_crossed(monitor):
    result = monitor.update(monitor_inputs(
        {"t": threshold(upper=10.0)}, {"Acetate": 2.0, "D-Glucose": 40.0, "E.coli": 0.5}))
    assert result["new_species"]["_add"] == {}
    assert result["new_species"]["_remove"] == []


def test_species_added_when_upper_threshold_crossed(monitor):
    result = monitor.update(monitor_inputs(
        {"t": threshold(upper=10.0)}, {"Acetate": 15.0, "D-Glucose": 40.0, "E.coli": 0.5}))

    added = result["new_species"]["_add"]
    assert set(added) == {"E.coli 2"}
    assert added["E.coli 2"]["address"] == "local:dFBA"
    assert added["E.coli 2"]["config"]["name"] == "E.coli 2"


def test_species_added_when_lower_threshold_crossed(monitor):
    """A lower-only threshold must fire.

    Regression test: unset `maybe` bounds are absent from the store rather than
    present as None, so indexing `range["lower"]` raised KeyError and no
    threshold could ever be evaluated.
    """
    result = monitor.update(monitor_inputs(
        {"t": threshold(lower=5.0)}, {"Acetate": 1.0, "D-Glucose": 40.0, "E.coli": 0.5}))
    assert set(result["new_species"]["_add"]) == {"E.coli 2"}


def test_lower_threshold_not_crossed(monitor):
    result = monitor.update(monitor_inputs(
        {"t": threshold(lower=5.0)}, {"Acetate": 8.0, "D-Glucose": 40.0, "E.coli": 0.5}))
    assert result["new_species"]["_add"] == {}


def test_new_species_mass_is_taken_from_the_parent(monitor):
    result = monitor.update(monitor_inputs(
        {"t": threshold(upper=10.0, mass=0.25)},
        {"Acetate": 15.0, "D-Glucose": 40.0, "E.coli": 0.5}))

    counts = result["counts"]
    assert counts["_add"] == {"E.coli 2": 0.25}
    assert counts["E.coli"] == -0.25, "parent must give up the mass it seeded"


def test_new_species_gets_a_results_store(monitor):
    result = monitor.update(monitor_inputs(
        {"t": threshold(upper=10.0)}, {"Acetate": 15.0, "D-Glucose": 40.0, "E.coli": 0.5}))
    added = result["dfba_results"]["_add"]["E.coli 2"]
    # substrates carried over from the parent, plus the new species' own biomass
    assert set(added) == {"D-Glucose", "Acetate", "E.coli 2"}


def test_existing_species_is_not_added_twice(monitor):
    result = monitor.update(monitor_inputs(
        {"t": threshold(upper=10.0)},
        {"Acetate": 15.0, "D-Glucose": 40.0, "E.coli": 0.5, "E.coli 2": 0.25},
        species=("E.coli", "E.coli 2")))
    assert result["new_species"]["_add"] == {}


def test_remove_threshold_marks_species_for_removal(monitor):
    result = monitor.update(monitor_inputs(
        {"t": threshold(upper=10.0, kind="remove")},
        {"Acetate": 15.0, "D-Glucose": 40.0, "E.coli": 0.5}))
    assert result["new_species"]["_remove"] == ["E.coli 2"]
    assert result["dfba_results"]["_remove"] == ["E.coli 2"]


def test_env_monitor_spec_wiring():
    spec = get_env_monitor_spec(interval=0.1)
    assert spec["_type"] == "step"
    assert spec["address"] == "local:EnvironmentMonitor"
    assert spec["inputs"]["thresholds"] == [THRESHOLDS]
    assert spec["outputs"]["counts"] == [SHARED_ENVIRONMENT, "counts"]


# --- integration -------------------------------------------------------------

def test_monitor_spawns_a_second_species_in_a_composite(core):
    """Acetate accumulates past the threshold, which spawns a daughter strain."""
    spec = build_spec(interval=0.1)
    spec[THRESHOLDS] = {"Acetate": threshold(upper=10.0)}
    spec["monitor"] = get_env_monitor_spec(interval=0.1)

    sim, results = run_spec(core, spec, 10)

    assert set(sim.state[SPECIES_STORE]) == {"E.coli", "E.coli 2"}
    assert "E.coli 2" in sim.state[DFBA_RESULTS]

    # the daughter only appears once seeded, so read the final timepoint
    final = results[-1]["shared_environment"]["concentrations"]
    assert final["E.coli 2"] > 0, "daughter strain should grow after being seeded"
    # it was seeded from the parent's biomass, not created from nothing
    assert final["Acetate"] > 10.0


def test_monitor_leaves_composite_alone_below_threshold(core):
    spec = build_spec(interval=0.1)
    spec[THRESHOLDS] = {"Acetate": threshold(upper=1e6)}
    spec["monitor"] = get_env_monitor_spec(interval=0.1)

    sim, _ = run_spec(core, spec, 5)
    assert set(sim.state[SPECIES_STORE]) == {"E.coli"}
