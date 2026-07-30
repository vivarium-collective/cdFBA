"""Tests for the environment-manipulating processes and steps.

Covers UpdateEnvironment, StaticConcentration, WaveFunction and Injector, both
as isolated units and wired into a running composite.
"""
from math import isclose, pi, sin

import pytest

from cdFBA.processes.dfba import (
    UpdateEnvironment, StaticConcentration, WaveFunction, Injector)
from cdFBA.utils import get_static_spec, get_wave_spec, get_injector_spec
from tests.conftest import build_spec, env, run_spec, series


def counts_update(result):
    return result["shared_environment"]["counts"]


# --- UpdateEnvironment -------------------------------------------------------

@pytest.fixture
def update_step(core):
    return UpdateEnvironment({}, core)


def test_update_environment_ports(update_step):
    assert set(update_step.inputs()) == {"shared_environment", "species_updates"}
    assert set(update_step.outputs()) == {"shared_environment"}


def test_update_environment_passes_through_a_single_species(update_step):
    shared = env({"D-Glucose": 10.0, "Acetate": 0.0}, volume=1.0)
    result = update_step.update({
        "shared_environment": shared,
        "species_updates": {"A": {"D-Glucose": -4.0, "Acetate": 2.0}},
    })
    assert counts_update(result) == {"D-Glucose": -4.0, "Acetate": 2.0}


def test_update_environment_sums_across_species(update_step):
    shared = env({"D-Glucose": 100.0, "Acetate": 0.0}, volume=1.0)
    result = update_step.update({
        "shared_environment": shared,
        "species_updates": {
            "A": {"D-Glucose": -4.0, "Acetate": 2.0},
            "B": {"D-Glucose": -3.0, "Acetate": 1.0},
        },
    })
    assert counts_update(result) == {"D-Glucose": -7.0, "Acetate": 3.0}


def test_update_environment_reports_zero_for_untouched_substrates(update_step):
    shared = env({"D-Glucose": 10.0, "Acetate": 5.0}, volume=1.0)
    result = update_step.update({
        "shared_environment": shared,
        "species_updates": {"A": {"D-Glucose": -1.0}},
    })
    assert counts_update(result)["Acetate"] == 0


def test_update_environment_clamps_a_single_species_at_zero(update_step):
    shared = env({"D-Glucose": 10.0, "Acetate": 0.0}, volume=1.0)
    result = update_step.update({
        "shared_environment": shared,
        "species_updates": {"A": {"D-Glucose": -15.0, "Acetate": 3.0}},
    })
    # glucose is drawn down to exactly zero, never below
    assert counts_update(result)["D-Glucose"] == -10.0


def test_update_environment_clamp_still_applies_other_substrates(update_step):
    """Hitting the clamp on one substrate must not drop the species' others.

    Regression test: the loop used to `break` out of the substrate loop as soon
    as one substrate was clamped, silently discarding every remaining exchange
    for that species.
    """
    shared = env({"D-Glucose": 10.0, "Acetate": 0.0}, volume=1.0)
    result = update_step.update({
        "shared_environment": shared,
        "species_updates": {"A": {"D-Glucose": -15.0, "Acetate": 3.0}},
    })
    assert counts_update(result)["Acetate"] == 3.0


@pytest.mark.parametrize("run", range(5))  # species order is randomised
def test_update_environment_clamps_across_multiple_species(update_step, run):
    """Combined demand above what is available must still clamp at zero.

    Regression test: each species used to be checked against the *original*
    amount rather than what was left, so two species could jointly drive a
    substrate negative.
    """
    shared = env({"D-Glucose": 10.0, "Acetate": 0.0}, volume=1.0)
    result = update_step.update({
        "shared_environment": shared,
        "species_updates": {
            "A": {"D-Glucose": -6.0, "Acetate": 3.0},
            "B": {"D-Glucose": -6.0, "Acetate": 3.0},
        },
    })
    update = counts_update(result)
    assert 10.0 + update["D-Glucose"] == pytest.approx(0.0)
    assert update["Acetate"] == 6.0


def test_update_environment_handles_no_species(update_step):
    shared = env({"D-Glucose": 10.0}, volume=1.0)
    result = update_step.update({"shared_environment": shared, "species_updates": {}})
    assert counts_update(result) == {"D-Glucose": 0}


# --- StaticConcentration -----------------------------------------------------

def test_static_concentration_returns_delta_to_target(core):
    process = StaticConcentration({"substrate_concentrations": {"D-Glucose": 40.0}}, core)
    shared = env({"D-Glucose": 30.0}, volume=2.0)  # concentration 15.0
    result = process.update({"shared_environment": shared, "global_time": 0.0}, 1.0)
    # target count is 40 * 2 = 80, current count is 30, so the delta is 50
    assert counts_update(result) == {"D-Glucose": 50.0}


def test_static_concentration_is_a_noop_when_already_at_target(core):
    process = StaticConcentration({"substrate_concentrations": {"D-Glucose": 15.0}}, core)
    shared = env({"D-Glucose": 30.0}, volume=2.0)
    result = process.update({"shared_environment": shared, "global_time": 0.0}, 1.0)
    assert counts_update(result)["D-Glucose"] == pytest.approx(0.0)


def test_static_concentration_holds_glucose_flat_in_a_composite(core):
    spec = build_spec()
    spec["StaticConcentration"] = get_static_spec(
        {"substrate_concentrations": {"D-Glucose": 40.0}}, interval=1.0)
    _, results = run_spec(core, spec, 10)

    glucose = series(results, "D-Glucose")
    assert all(g == pytest.approx(40.0) for g in glucose[1:])
    # the culture still grows, fed by the held-constant glucose
    biomass = series(results, "E.coli")
    assert biomass[-1] > biomass[0]


# --- WaveFunction ------------------------------------------------------------

def test_wave_function_tracks_the_target_wave(core):
    params = {"D-Glucose": {"amplitude": 5.0, "angular_frequency": 1.0,
                            "base_concentration": 10.0, "phase_shift": 0.0}}
    process = WaveFunction({"substrate_params": params}, core)
    shared = env({"D-Glucose": 20.0}, volume=2.0)  # current count 20

    t = pi / 2  # sin = 1, so target concentration is 15
    result = process.update({"shared_environment": shared, "global_time": t}, 1.0)
    assert counts_update(result)["D-Glucose"] == pytest.approx(15.0 * 2 - 20.0)


def test_wave_function_at_zero_phase_targets_base_concentration(core):
    params = {"D-Glucose": {"amplitude": 5.0, "angular_frequency": 1.0,
                            "base_concentration": 10.0, "phase_shift": 0.0}}
    process = WaveFunction({"substrate_params": params}, core)
    shared = env({"D-Glucose": 20.0}, volume=2.0)
    result = process.update({"shared_environment": shared, "global_time": 0.0}, 1.0)
    # target is 10 * 2 = 20, which is exactly the current count
    assert counts_update(result)["D-Glucose"] == pytest.approx(0.0)


def test_wave_function_does_not_drive_negative_targets(core):
    """When the wave dips to or below zero the process makes no change."""
    params = {"D-Glucose": {"amplitude": 10.0, "angular_frequency": 1.0,
                            "base_concentration": 1.0, "phase_shift": 0.0}}
    process = WaveFunction({"substrate_params": params}, core)
    shared = env({"D-Glucose": 20.0}, volume=2.0)
    t = 3 * pi / 2  # sin = -1, target concentration is 1 - 10 = -9
    result = process.update({"shared_environment": shared, "global_time": t}, 1.0)
    assert counts_update(result)["D-Glucose"] == 0


def test_wave_function_oscillates_in_a_composite(core):
    spec = build_spec(glucose=10.0)
    spec["WaveFunction"] = get_wave_spec({"substrate_params": {
        "D-Glucose": {"amplitude": 5.0, "angular_frequency": pi / 5,
                      "base_concentration": 10.0, "phase_shift": 0.0}}},
        interval=1.0)
    _, results = run_spec(core, spec, 20)

    glucose = series(results, "D-Glucose")
    # a full period is 10 steps, so the series must both rise and fall
    assert max(glucose) > min(glucose) + 1.0
    assert max(glucose) <= 15.0 + 1e-6


# --- Injector ----------------------------------------------------------------

def make_injector(core, params):
    return Injector({"injection_params": params}, core)


@pytest.mark.parametrize("time,expected", [
    (0.0, None),   # nothing is injected at t=0
    (2.5, None),
    (5.0, 80.0),
    (7.5, None),
    (10.0, 80.0),
    (15.0, 80.0),
])
def test_injector_fires_on_a_fixed_interval(core, time, expected):
    process = make_injector(core, {"D-Glucose": {"amount": 80.0, "interval": 5.0}})
    shared = env({"D-Glucose": 10.0}, volume=1.0)
    update = counts_update(
        process.update({"shared_environment": shared, "global_time": time}, 1.0))
    assert update.get("D-Glucose") == expected


@pytest.mark.parametrize("time,expected", [
    (2.0, 10.0),   # first interval
    (5.0, 10.0),   # 2 + 3
    (7.0, 10.0),   # 2 + 3 + 2
    (10.0, 10.0),  # 2 + 3 + 2 + 3
    (3.0, None),
    (4.0, None),
    (6.0, None),
])
def test_injector_cycles_through_a_list_of_intervals(core, time, expected):
    process = make_injector(core, {"D-Glucose": {"amount": 10.0, "interval": [2.0, 3.0]}})
    shared = env({"D-Glucose": 0.0}, volume=1.0)
    update = counts_update(
        process.update({"shared_environment": shared, "global_time": time}, 1.0))
    assert update.get("D-Glucose") == expected


def test_injector_handles_multiple_substrates(core):
    process = make_injector(core, {
        "D-Glucose": {"amount": 80.0, "interval": 5.0},
        "Acetate": {"amount": 5.0, "interval": 2.0},
    })
    shared = env({"D-Glucose": 0.0, "Acetate": 0.0}, volume=1.0)
    at_ten = counts_update(
        process.update({"shared_environment": shared, "global_time": 10.0}, 1.0))
    assert at_ten == {"D-Glucose": 80.0, "Acetate": 5.0}

    at_four = counts_update(
        process.update({"shared_environment": shared, "global_time": 4.0}, 1.0))
    assert at_four == {"Acetate": 5.0}


def test_injector_replenishes_glucose_in_a_composite(core):
    spec = build_spec()
    spec["Injector"] = get_injector_spec({"injection_params": {
        "D-Glucose": {"amount": 80.0, "interval": 5.0}}}, interval=1.0)
    _, results = run_spec(core, spec, 20)

    glucose = series(results, "D-Glucose")
    plain = series(run_spec(core, build_spec(), 20)[1], "D-Glucose")

    # injections keep more glucose around than an un-fed run
    assert sum(glucose) > sum(plain)
    assert series(results, "E.coli")[-1] > series(results, "E.coli")[0]
