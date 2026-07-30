"""Whole-composite integration tests, including co-culture runs.

The `bigg` marked test pulls models from the BiGG database and needs network
access; deselect it with ``-m "not bigg"``.
"""
import pytest
from process_bigraph import Composite
from process_bigraph.emitter import gather_emitter_results, emitter_from_wires

from cdFBA.utils import SHARED_ENVIRONMENT, SPECIES_STORE, DFBA_RESULTS
from cdFBA.utils import make_cdfba_composite, set_concentration, set_kinetics
from tests.conftest import build_spec, run_spec, series, EXCHANGES, KINETICS, BOUNDS


def test_two_species_share_one_environment(core):
    spec = build_spec(species=("E.coli", "S.other"))
    sim, results = run_spec(core, spec, 10)

    assert set(sim.state[SPECIES_STORE]) == {"E.coli", "S.other"}

    glucose = series(results, "D-Glucose")
    assert glucose[0] > glucose[-1], "the co-culture should draw down glucose"
    assert all(g >= -1e-9 for g in glucose), "glucose must never go negative"

    for name in ("E.coli", "S.other"):
        biomass = series(results, name)
        assert biomass[-1] > biomass[0], f"{name} should grow"


def test_two_species_deplete_glucose_faster_than_one(core):
    def exhausted_at(results):
        glucose = series(results, "D-Glucose")
        return next(i for i, g in enumerate(glucose) if g <= 1e-9)

    _, solo = run_spec(core, build_spec(species=("E.coli",)), 10)
    _, duo = run_spec(core, build_spec(species=("E.coli", "S.other")), 10)

    # both cultures run out of glucose, but two consumers get there sooner
    assert exhausted_at(duo) < exhausted_at(solo)


def test_glucose_is_never_driven_negative_when_scarce(core):
    """With little glucose and two consumers, the clamp must hold at zero."""
    spec = build_spec(species=("E.coli", "S.other"), glucose=1.0)
    _, results = run_spec(core, spec, 10)

    glucose = series(results, "D-Glucose")
    assert min(glucose) >= -1e-9
    assert glucose[-1] == pytest.approx(0.0, abs=1e-6), "glucose should be exhausted"


def test_counts_and_concentrations_stay_consistent(core):
    """The volumetric type must keep counts == concentrations * volume."""
    spec = build_spec(volume=2.0)
    _, results = run_spec(core, spec, 5)

    sim = Composite({"state": build_spec(volume=2.0)}, core=core)
    sim.run(5)
    store = sim.state[SHARED_ENVIRONMENT]
    volume = store["volume"]
    for substrate, count in store["counts"].items():
        assert store["concentrations"][substrate] == pytest.approx(count / volume)


def test_emitter_records_one_result_per_step(core, spec):
    _, results = run_spec(core, spec, 20)
    assert len(results) == 21
    assert [r["global_time"] for r in results] == [float(t) for t in range(21)]


@pytest.mark.bigg
def test_bigg_cocultures_run(core):
    """End-to-end run on the genome-scale BiGG models used by the original demo."""
    model_dict = {"E.coli": "iAF1260", "S.flexneri": "iSFxv_1172"}
    spec = make_cdfba_composite(
        model_dict, medium_type=None, exchanges=EXCHANGES, volume=2, interval=1.0)

    for name in model_dict:
        spec[SPECIES_STORE][name]["config"]["bounds"] = dict(BOUNDS)
        set_kinetics(name, spec, KINETICS)
    set_concentration(spec, {"Acetate": 0, "D-Glucose": 40})
    spec["emitter"] = emitter_from_wires({
        "global_time": ["global_time"],
        "shared_environment": [SHARED_ENVIRONMENT],
    })

    sim = Composite({"state": spec}, core=core)
    sim.run(20)
    results = gather_emitter_results(sim)[("emitter",)]

    assert len(results) == 21
    assert series(results, "D-Glucose")[2] > series(results, "D-Glucose")[4]
    assert series(results, "E.coli")[4] > series(results, "E.coli")[2]
