"""Shared fixtures and helpers for the cdFBA test suite.

Tests default to cobra's bundled ``textbook`` (e_coli_core) model so the suite
runs fast and without network access. Tests that pull models from BiGG are
marked ``bigg`` and can be deselected with ``-m "not bigg"``.
"""
import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results

from cdFBA.data_types import register_types
from cdFBA.utils import SHARED_ENVIRONMENT, SPECIES_STORE
from cdFBA.utils import make_cdfba_composite, set_concentration, set_kinetics

TEXTBOOK = "textbook"
EXCHANGES = ["EX_glc__D_e", "EX_ac_e"]
KINETICS = {"D-Glucose": (0.02, 15), "Acetate": (0.5, 7)}
BOUNDS = {
    "EX_o2_e": {"lower": -2, "upper": None},
    "ATPM": {"lower": 1, "upper": 1},
}


@pytest.fixture(scope="session")
def core():
    """A core with the cdFBA types and processes registered."""
    return register_types(allocate_core())


def build_spec(species=("E.coli",), volume=2.0, interval=1.0,
               glucose=40.0, acetate=0.0, emitter=True):
    """Build a textbook-model cdFBA composite spec with an emitter attached."""
    spec = make_cdfba_composite(
        {name: TEXTBOOK for name in species},
        medium_type=None,
        exchanges=EXCHANGES,
        volume=volume,
        interval=interval,
    )
    for name in species:
        spec[SPECIES_STORE][name]["config"]["bounds"] = dict(BOUNDS)
        set_kinetics(name, spec, KINETICS)
    set_concentration(spec, {"D-Glucose": glucose, "Acetate": acetate})
    if emitter:
        spec["emitter"] = emitter_from_wires({
            "global_time": ["global_time"],
            "shared_environment": [SHARED_ENVIRONMENT],
        })
    return spec


@pytest.fixture
def spec():
    return build_spec()


def run_spec(core, spec, duration):
    """Run a spec as a Composite and return ``(sim, emitted_results)``."""
    sim = Composite({"state": spec}, core=core)
    sim.run(duration)
    return sim, gather_emitter_results(sim)[("emitter",)]


def series(results, substrate):
    """Extract the concentration time-series for one substrate."""
    return [r["shared_environment"]["concentrations"][substrate] for r in results]


def env(counts, volume=2.0):
    """Build a volumetric store value from a counts dict."""
    return {
        "volume": volume,
        "counts": dict(counts),
        "concentrations": {k: v / volume for k, v in counts.items()},
    }
