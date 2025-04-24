from pprint import pprint
import numpy as np

from process_bigraph import ProcessTypes, Process, Step, Composite
from process_bigraph.emitter import gather_emitter_results

from cdFBA.utils import SHARED_ENVIRONMENT, SPECIES_STORE, THRESHOLDS, DFBA_RESULTS, FIELDS
from cdFBA.utils import get_single_dfba_spec, set_concentration, make_cdfba_composite, set_kinetics
from cdFBA.processes.dfba import dFBA, UpdateEnvironment

from matplotlib import pyplot as plt

def create_spatial(dims, distance):
    """Creates a spec for shared environments in Euclidean Space

    Parameters:
        dims: list of int, number of compartments in each spatial dimension [x, y, z]
        distance: float, distance between neighboring voxels
    """

    fields = {
        FIELDS: {}
    }
    field = 0
    for i in np.arange(distance/2, distance*dims[0], distance):
        for j in np.arange(distance/2, distance*dims[1], distance):
            if dims[2] != 0:
                for k in np.arange(distance/2, distance*dims[2], distance):
                    fields[FIELDS][f"[{field}]"] = {}
                    fields[FIELDS][f"[{field}]"]["location"] = [i, j, k]
                    field += 1
            else:
                fields[FIELDS][f"[{field}]"] = {}
                k=0
                fields[FIELDS][f"[{field}]"]["location"] = [i, j, k]
                field += 1

    return fields

if __name__ == "__main__":
    fields = create_spatial(dims=[2, 2, 2], distance=1)
    pprint(fields)