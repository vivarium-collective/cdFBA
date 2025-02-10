from cdFBA.processes.dfba import DFBA, UpdateEnvironment, Chemostat, Injector, WaveFunction
from cdFBA.utils import DFBAconfig, model_from_file, get_objective_reaction, get_injector_spec, get_wave_spec, get_chemo_spec
from cdFBA.utils import dfba_config_from_model, get_single_dfba_spec, dfba_config, environment_spec, initial_environment
import cobra

from process_bigraph.composite import ProcessTypes
from process_bigraph import Composite

from cdFBA import register_types

from matplotlib import pyplot as plt

from vivarium import Vivarium

processes = {
    'DFBA': DFBA,
    'Chemostat': Chemostat,
    'Injector': Injector,
    'WaveFunction': WaveFunction,
    'UpdateEnvironment': UpdateEnvironment,
}

v = Vivarium(processes=processes)
print(v.process_interface('DFBA'))