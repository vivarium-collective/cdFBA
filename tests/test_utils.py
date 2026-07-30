"""Tests for the model-introspection and spec-building helpers in cdFBA.utils."""
import pytest

from cdFBA.utils import SHARED_ENVIRONMENT, SPECIES_STORE, DFBA_RESULTS
from cdFBA.utils import (
    model_from_file, get_model_dict, get_exchanges, get_substrates,
    get_reaction_map, get_kinetics, get_bounds, get_objective_reaction,
    dfba_config, get_single_dfba_spec, make_cdfba_composite,
    get_combined_exchanges, get_initial_counts, initial_environment,
    set_counts, set_concentration, set_kinetics,
    environment_spec, get_static_spec, get_wave_spec, get_injector_spec,
)
from tests.conftest import TEXTBOOK, EXCHANGES


@pytest.fixture(scope="module")
def model():
    return model_from_file(TEXTBOOK)


def test_model_from_file_loads_bundled_model(model):
    assert model.id == "e_coli_core"


def test_model_from_file_rejects_non_string():
    with pytest.raises(ValueError):
        model_from_file(42)


def test_get_model_dict_loads_each_entry():
    models = get_model_dict({"a": TEXTBOOK, "b": TEXTBOOK})
    assert set(models) == {"a", "b"}
    assert all(m.id == "e_coli_core" for m in models.values())


def test_get_exchanges_medium_types(model):
    default = get_exchanges(model, medium_type="default")
    every = get_exchanges(model, medium_type="exchange")
    minimal = get_exchanges(model, medium_type="minimal")

    # "exchange" is the superset; the default medium is a subset of it
    assert set(default).issubset(set(every))
    assert set(minimal).issubset(set(every))
    assert len(every) > len(default)


def test_get_exchanges_rejects_unknown_medium_type(model):
    with pytest.raises(ValueError):
        get_exchanges(model, medium_type="nonsense")


def test_get_substrates_and_reaction_map_agree(model):
    substrates = get_substrates(model, EXCHANGES)
    reaction_map = get_reaction_map(model, EXCHANGES)

    assert substrates == ["D-Glucose", "Acetate"]
    assert reaction_map == {"D-Glucose": "EX_glc__D_e", "Acetate": "EX_ac_e"}
    assert set(reaction_map) == set(substrates)


def test_reaction_map_skips_reactions_absent_from_model(model):
    reaction_map = get_reaction_map(model, EXCHANGES + ["EX_not_a_reaction_e"])
    assert set(reaction_map) == {"D-Glucose", "Acetate"}


def test_get_kinetics_keys_match_substrates(model):
    kinetics = get_kinetics(model, EXCHANGES)
    assert set(kinetics) == set(get_substrates(model, EXCHANGES))
    assert all(len(params) == 2 for params in kinetics.values())


def test_get_bounds_is_keyed_by_reaction_id(model):
    bounds = get_bounds(get_reaction_map(model, EXCHANGES), upper=10, lower=-10)
    assert bounds == {
        "EX_glc__D_e": {"lower": -10, "upper": 10},
        "EX_ac_e": {"lower": -10, "upper": 10},
    }


def test_get_objective_reaction(model):
    assert get_objective_reaction(model) == "Biomass_Ecoli_core"


def test_dfba_config_defaults_are_self_consistent():
    """Default kinetics must cover every substrate in the default reaction map.

    Regression test: the defaults previously used lowercase substrate names
    ("glucose") that never matched the reaction map ("D-Glucose"), so a config
    built purely from defaults raised KeyError on the first update.
    """
    config = dfba_config(model_file=TEXTBOOK)
    assert set(config["reaction_map"]).issubset(set(config["kinetics"]))


def test_dfba_config_defaults_name_to_model_id():
    assert dfba_config(model_file=TEXTBOOK)["name"] == "e_coli_core"


def test_dfba_config_passes_through_overrides():
    config = dfba_config(
        model_file=TEXTBOOK,
        name="strain",
        kinetics={"D-Glucose": (1.0, 2.0)},
        reaction_map={"D-Glucose": "EX_glc__D_e"},
        bounds={},
        medium={"EX_glc__D_e": 10},
    )
    assert config["name"] == "strain"
    assert config["bounds"] == {}
    assert config["medium"] == {"EX_glc__D_e": 10}


def test_get_single_dfba_spec_wiring():
    spec = get_single_dfba_spec(model_file=TEXTBOOK, name="E.coli", interval=0.5)

    assert spec["_type"] == "process"
    assert spec["address"] == "local:dFBA"
    assert spec["interval"] == 0.5
    assert spec["inputs"]["shared_environment"] == ["..", SHARED_ENVIRONMENT]
    assert spec["outputs"]["dfba_update"] == ["..", DFBA_RESULTS, "E.coli"]


def test_initial_environment_derives_concentrations():
    store = initial_environment(volume=4, initial_counts={"D-Glucose": 40, "E.coli": 2})
    assert store["volume"] == 4
    assert store["concentrations"] == {"D-Glucose": 10, "E.coli": 0.5}


def test_initial_environment_requires_counts_or_species():
    with pytest.raises(ValueError):
        initial_environment(volume=1)


def test_get_initial_counts_covers_substrates_and_biomass():
    counts = get_initial_counts(
        {"E.coli": TEXTBOOK}, biomass=0.5, initial_value=20, exchanges=EXCHANGES)
    assert counts == {"D-Glucose": 20, "Acetate": 20, "E.coli": 0.5}


def test_get_initial_counts_requires_exchanges():
    with pytest.raises(ValueError):
        get_initial_counts({"E.coli": TEXTBOOK})


def test_get_combined_exchanges_is_a_union(model):
    combined = get_combined_exchanges({"a": model, "b": model}, medium_type="default")
    assert set(combined) == set(get_exchanges(model, medium_type="default"))


# --- spec mutation helpers ---------------------------------------------------

def test_set_concentration_updates_counts(spec):
    set_concentration(spec, {"D-Glucose": 10.0})
    assert spec[SHARED_ENVIRONMENT]["concentrations"]["D-Glucose"] == 10.0
    # volume is 2.0, so counts must be twice the concentration
    assert spec[SHARED_ENVIRONMENT]["counts"]["D-Glucose"] == 20.0


def test_set_counts_updates_concentrations(spec):
    set_counts(spec, {"D-Glucose": 30.0})
    assert spec[SHARED_ENVIRONMENT]["counts"]["D-Glucose"] == 30.0
    assert spec[SHARED_ENVIRONMENT]["concentrations"]["D-Glucose"] == 15.0


@pytest.mark.parametrize("mutate", [set_counts, set_concentration])
def test_setters_reject_unknown_substrate(spec, mutate):
    with pytest.raises(ValueError, match="not in shared environment"):
        mutate(spec, {"Nonexistent": 1.0})


def test_set_kinetics_writes_into_species_config(spec):
    set_kinetics("E.coli", spec, {"D-Glucose": (1.5, 9.0)})
    kinetics = spec[SPECIES_STORE]["E.coli"]["config"]["kinetics"]
    assert kinetics["D-Glucose"] == (1.5, 9.0)


def test_set_kinetics_rejects_unknown_substrate(spec):
    with pytest.raises(ValueError, match="not in shared environment"):
        set_kinetics("E.coli", spec, {"Nonexistent": (1.0, 1.0)})


# --- composite assembly ------------------------------------------------------

def test_make_cdfba_composite_structure():
    spec = make_cdfba_composite(
        {"E.coli": TEXTBOOK, "S.other": TEXTBOOK},
        medium_type=None, exchanges=EXCHANGES, volume=2, interval=1.0)

    assert set(spec[SPECIES_STORE]) == {"E.coli", "S.other"}
    assert set(spec[DFBA_RESULTS]) == {"E.coli", "S.other"}
    assert spec["update environment"]["address"] == "local:UpdateEnvironment"
    assert spec[SHARED_ENVIRONMENT]["volume"] == 2

    # every species biomass and shared substrate is present in the environment
    counts = spec[SHARED_ENVIRONMENT]["counts"]
    assert {"D-Glucose", "Acetate", "E.coli", "S.other"} == set(counts)

    # each species' results store is seeded with its substrates plus its biomass
    assert set(spec[DFBA_RESULTS]["E.coli"]) == {"D-Glucose", "Acetate", "E.coli"}


def test_make_cdfba_composite_requires_medium_type_or_exchanges():
    with pytest.raises(ValueError, match="Must provide"):
        make_cdfba_composite({"E.coli": TEXTBOOK}, medium_type=None, exchanges=None)


def test_make_cdfba_composite_rejects_both_medium_type_and_exchanges():
    with pytest.raises(ValueError, match="only on"):
        make_cdfba_composite(
            {"E.coli": TEXTBOOK}, medium_type="default", exchanges=EXCHANGES)


def test_make_cdfba_composite_from_medium_type():
    spec = make_cdfba_composite(
        {"E.coli": TEXTBOOK}, medium_type="default", exchanges=None, volume=1)
    reaction_map = spec[SPECIES_STORE]["E.coli"]["config"]["reaction_map"]
    assert "D-Glucose" in reaction_map
    assert set(reaction_map).issubset(set(spec[SHARED_ENVIRONMENT]["counts"]))


# --- environmental process specs --------------------------------------------

def test_environment_spec_wiring():
    step = environment_spec()
    assert step["_type"] == "step"
    assert step["inputs"]["species_updates"] == [DFBA_RESULTS]
    assert step["outputs"]["shared_environment"] == [SHARED_ENVIRONMENT]


@pytest.mark.parametrize("builder,address", [
    (get_static_spec, "local:StaticConcentration"),
    (get_wave_spec, "local:WaveFunction"),
    (get_injector_spec, "local:Injector"),
])
def test_environmental_specs_wire_time_and_environment(builder, address):
    spec = builder({"anything": {}}, interval=0.25)
    assert spec["address"] == address
    assert spec["interval"] == 0.25
    assert spec["inputs"]["global_time"] == ["global_time"]
    assert spec["inputs"]["shared_environment"] == [SHARED_ENVIRONMENT]
    assert spec["outputs"]["shared_environment"] == [SHARED_ENVIRONMENT]


@pytest.mark.parametrize("builder", [get_static_spec, get_wave_spec, get_injector_spec])
def test_environmental_specs_require_config(builder):
    with pytest.raises(ValueError, match="provide config"):
        builder()
