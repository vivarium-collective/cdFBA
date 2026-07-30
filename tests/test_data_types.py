"""Tests for the custom cdFBA schema types, especially `volumetric`."""
import pytest


BASE = {
    "volume": 2.0,
    "counts": {"A": 10.0, "B": 4.0},
    "concentrations": {"A": 5.0, "B": 2.0},
}


def apply_update(core, update, base=None):
    """Apply an update to a fresh `volumetric` value and return the result."""
    schema, state, _ = core.realize("volumetric", base or BASE)
    applied = core.apply(schema, state, update)
    return applied[0] if isinstance(applied, tuple) else applied


def test_registered_types_are_available(core):
    for type_name in ("volumetric", "bounds", "threshold", "dfba_changes"):
        assert core.access(type_name) is not None


def test_realize_keeps_all_volumetric_fields(core):
    _, state, _ = core.realize("volumetric", BASE)
    assert state["volume"] == 2.0
    assert state["counts"] == {"A": 10.0, "B": 4.0}
    assert state["concentrations"] == {"A": 5.0, "B": 2.0}


def test_counts_update_recomputes_concentrations(core):
    result = apply_update(core, {"counts": {"A": 5.0}})
    # counts apply additively, concentrations are derived from counts / volume
    assert result["counts"] == {"A": 15.0, "B": 4.0}
    assert result["concentrations"] == {"A": 7.5, "B": 2.0}
    assert result["volume"] == 2.0


def test_volume_update_recomputes_concentrations(core):
    result = apply_update(core, {"volume": 2.0})
    assert result["volume"] == 4.0
    assert result["counts"] == {"A": 10.0, "B": 4.0}
    assert result["concentrations"] == {"A": 2.5, "B": 1.0}


def test_concentrations_update_recomputes_counts(core):
    """A concentrations update must take effect and back-propagate to counts.

    Regression test: the dispatch guard checked for the singular key
    "concentration", which never matched the "concentrations" field, so every
    concentration update was silently swallowed and the state left unchanged.
    """
    result = apply_update(core, {"concentrations": {"A": 9.0}})
    assert result["concentrations"]["A"] == 14.0
    assert result["counts"]["A"] == 28.0
    assert result["volume"] == 2.0


def test_concentrations_and_volume_update_together_is_rejected(core):
    with pytest.raises(ValueError, match="volume and concentration"):
        apply_update(core, {"concentrations": {"A": 1.0}, "volume": 1.0})


def test_concentrations_and_counts_update_together_is_rejected(core):
    with pytest.raises(ValueError, match="count and concentration"):
        apply_update(core, {"concentrations": {"A": 1.0}, "counts": {"A": 1.0}})


def test_empty_update_is_a_noop(core):
    result = apply_update(core, {})
    assert result["counts"] == BASE["counts"]
    assert result["concentrations"] == BASE["concentrations"]
    assert result["volume"] == BASE["volume"]


def test_unset_maybe_bounds_are_absent_rather_than_none(core):
    """`maybe` fields left unset drop out of the realized state entirely.

    Process code must therefore read optional fields with `.get()` rather than
    indexing and testing for None.
    """
    _, state, _ = core.realize("bounds", {"lower": None, "upper": 10.0})
    assert state["upper"] == 10.0
    assert "lower" not in state
