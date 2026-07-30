"""Tests for the dFBA process itself: construction, config handling, and updates."""
import pytest

from cdFBA.processes.dfba import dFBA
from cdFBA.utils import SHARED_ENVIRONMENT, dfba_config
from tests.conftest import TEXTBOOK, EXCHANGES, KINETICS, BOUNDS, env, run_spec, series

REACTION_MAP = {"D-Glucose": "EX_glc__D_e", "Acetate": "EX_ac_e"}


def make_config(**overrides):
    config = dfba_config(
        model_file=TEXTBOOK,
        name="E.coli",
        kinetics=dict(KINETICS),
        reaction_map=dict(REACTION_MAP),
        bounds=dict(BOUNDS),
    )
    config.update(overrides)
    return config


def make_process(core, **overrides):
    return dFBA(make_config(**overrides), core)


ENVIRONMENT = env({"D-Glucose": 80.0, "Acetate": 0.0, "E.coli": 1.0}, volume=2.0)


# --- construction / config ---------------------------------------------------

def test_biomass_identifier_is_discovered(core):
    assert make_process(core).biomass_identifier == "Biomass_Ecoli_core"


def test_bounds_from_config_are_applied_to_the_model(core):
    process = make_process(core)
    assert process.model.reactions.ATPM.bounds == (1, 1)
    # "upper": None must leave the model's own upper bound untouched
    assert process.model.reactions.EX_o2_e.lower_bound == -2
    assert process.model.reactions.EX_o2_e.upper_bound == 1000.0


@pytest.mark.parametrize("field", ["medium", "bounds", "changes"])
def test_optional_config_fields_may_be_absent(core, field):
    """Unset `maybe` config fields drop out of config entirely.

    Regression test: dFBA.__init__ indexed these fields directly, so an omitted
    (or None) medium/bounds/changes raised KeyError at construction.
    """
    config = make_config()
    del config[field]
    process = dFBA(config, core)
    assert process.model is not None


@pytest.mark.parametrize("field", ["medium", "bounds", "changes"])
def test_optional_config_fields_may_be_none(core, field):
    process = dFBA(make_config(**{field: None}), core)
    assert process.model is not None


def test_medium_is_applied_when_provided(core):
    medium = {"EX_glc__D_e": 5.0, "EX_o2_e": 1000.0, "EX_h_e": 1000.0,
              "EX_h2o_e": 1000.0, "EX_nh4_e": 1000.0, "EX_pi_e": 1000.0,
              "EX_co2_e": 1000.0}
    process = make_process(core, medium=medium)
    assert process.model.medium["EX_glc__D_e"] == 5.0


def test_changes_apply_reaction_knockout(core):
    process = make_process(core, changes={
        "gene_knockout": [], "reaction_knockout": ["PGI"],
        "bounds": {}, "kinetics": {}})
    assert process.model.reactions.PGI.bounds == (0, 0)


def test_changes_apply_gene_knockout(core):
    process = make_process(core, changes={
        "gene_knockout": ["b4025"], "reaction_knockout": [],
        "bounds": {}, "kinetics": {}})
    # b4025 is the only gene for PGI, so knocking it out closes the reaction
    assert process.model.reactions.PGI.bounds == (0, 0)


def test_changes_apply_bounds(core):
    process = make_process(core, changes={
        "gene_knockout": [], "reaction_knockout": [],
        "bounds": {"PGI": {"lower": -5.0, "upper": 5.0}}, "kinetics": {}})
    assert process.model.reactions.PGI.bounds == (-5.0, 5.0)


def test_changes_kinetics_override_config_kinetics(core):
    process = make_process(core, changes={
        "gene_knockout": [], "reaction_knockout": [],
        "bounds": {}, "kinetics": {"D-Glucose": (1.0, 99.0)}})
    assert process.config["kinetics"]["D-Glucose"] == (1.0, 99.0)
    # untouched substrates keep their original parameters
    assert process.config["kinetics"]["Acetate"] == KINETICS["Acetate"]


# --- ports -------------------------------------------------------------------

def test_ports(core):
    process = make_process(core)
    assert set(process.inputs()) == {"shared_environment", "current_update"}
    assert set(process.outputs()) == {"dfba_update"}


# --- update ------------------------------------------------------------------

def test_update_returns_a_delta_for_every_substrate_and_biomass(core):
    update = make_process(core).update(
        {"shared_environment": ENVIRONMENT, "current_update": {}}, 1.0)["dfba_update"]
    assert set(update) == {"D-Glucose", "Acetate", "E.coli"}


def test_update_grows_biomass_and_consumes_glucose(core):
    update = make_process(core).update(
        {"shared_environment": ENVIRONMENT, "current_update": {}}, 1.0)["dfba_update"]
    assert update["E.coli"] > 0, "biomass should increase on glucose"
    assert update["D-Glucose"] < 0, "glucose should be consumed"


def test_update_scales_with_interval(core):
    process = make_process(core)
    inputs = {"shared_environment": ENVIRONMENT, "current_update": {}}
    short = process.update(inputs, 1.0)["dfba_update"]
    long = process.update(inputs, 2.0)["dfba_update"]
    assert long["E.coli"] == pytest.approx(2 * short["E.coli"])
    assert long["D-Glucose"] == pytest.approx(2 * short["D-Glucose"])


def test_update_scales_with_current_biomass(core):
    process = make_process(core)
    single = process.update({"shared_environment": ENVIRONMENT, "current_update": {}}, 1.0)
    doubled_env = env({"D-Glucose": 80.0, "Acetate": 0.0, "E.coli": 2.0}, volume=2.0)
    doubled = process.update({"shared_environment": doubled_env, "current_update": {}}, 1.0)
    assert doubled["dfba_update"]["E.coli"] == pytest.approx(
        2 * single["dfba_update"]["E.coli"])


def test_uptake_is_michaelis_menten_limited(core):
    """Uptake flux is bounded by Vmax, and rises with substrate concentration."""
    process = make_process(core)
    km, vmax = KINETICS["D-Glucose"]

    scarce = env({"D-Glucose": 0.02, "Acetate": 0.0, "E.coli": 1.0}, volume=2.0)
    plenty = env({"D-Glucose": 200.0, "Acetate": 0.0, "E.coli": 1.0}, volume=2.0)

    low = process.update({"shared_environment": scarce, "current_update": {}}, 1.0)
    high = process.update({"shared_environment": plenty, "current_update": {}}, 1.0)

    assert abs(low["dfba_update"]["D-Glucose"]) < abs(high["dfba_update"]["D-Glucose"])
    assert abs(high["dfba_update"]["D-Glucose"]) <= vmax + 1e-6


def test_no_growth_without_glucose(core):
    starved = env({"D-Glucose": 0.0, "Acetate": 0.0, "E.coli": 1.0}, volume=2.0)
    update = make_process(core).update(
        {"shared_environment": starved, "current_update": {}}, 1.0)["dfba_update"]
    assert update["E.coli"] == pytest.approx(0.0, abs=1e-9)


# --- integration -------------------------------------------------------------

def test_composite_grows_and_depletes_glucose(core, spec):
    _, results = run_spec(core, spec, 10)

    assert len(results) == 11
    glucose = series(results, "D-Glucose")
    biomass = series(results, "E.coli")

    assert glucose[0] > glucose[-1], "glucose should be consumed over time"
    assert biomass[-1] > biomass[0], "biomass should grow over time"
    assert all(g >= -1e-9 for g in glucose), "glucose must never go negative"


def test_composite_conserves_biomass_monotonicity(core, spec):
    _, results = run_spec(core, spec, 10)
    biomass = series(results, "E.coli")
    assert all(b2 >= b1 - 1e-9 for b1, b2 in zip(biomass, biomass[1:]))
