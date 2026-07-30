"""Tests for the spatial field layout helper."""
import pytest

from cdFBA.processes.spatial import create_spatial
from cdFBA.utils import FIELDS


def locations(fields):
    return [tuple(float(c) for c in v["location"]) for v in fields[FIELDS].values()]


def test_two_dimensional_grid_is_flat():
    fields = create_spatial(dims=[2, 2, 0], distance=1.0)

    assert len(fields[FIELDS]) == 4
    assert locations(fields) == [
        (0.5, 0.5, 0.0), (0.5, 1.5, 0.0), (1.5, 0.5, 0.0), (1.5, 1.5, 0.0)]


def test_three_dimensional_grid_fills_the_volume():
    fields = create_spatial(dims=[2, 2, 2], distance=1.0)

    assert len(fields[FIELDS]) == 8
    assert (0.5, 0.5, 0.5) in locations(fields)
    assert (1.5, 1.5, 1.5) in locations(fields)


@pytest.mark.parametrize("dims,expected", [
    ([1, 1, 0], 1),
    ([3, 2, 0], 6),
    ([2, 3, 4], 24),
])
def test_field_count_matches_dimensions(dims, expected):
    assert len(create_spatial(dims=dims, distance=1.0)[FIELDS]) == expected


def test_voxels_are_centred_and_spaced_by_distance():
    fields = create_spatial(dims=[2, 1, 0], distance=4.0)
    xs = sorted(loc[0] for loc in locations(fields))
    # voxel centres sit half a step in, one full step apart
    assert xs == [2.0, 6.0]


def test_field_keys_are_unique_and_sequential():
    fields = create_spatial(dims=[2, 2, 0], distance=1.0)
    assert list(fields[FIELDS]) == ["[0]", "[1]", "[2]", "[3]"]
