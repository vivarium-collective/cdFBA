from cdFBA.processes.dfba import dFBA, UpdateEnvironment, StaticConcentration, Injector, WaveFunction
from cdFBA.processes.dfbalauncher import EnvironmentMonitor

def register_processes(core):
    core.register_link('dFBA', dFBA)
    core.register_link('UpdateEnvironment', UpdateEnvironment)
    core.register_link('StaticConcentration', StaticConcentration)
    core.register_link('Injector', Injector)
    core.register_link('WaveFunction', WaveFunction)
    core.register_link('EnvironmentMonitor', EnvironmentMonitor)
    
    return core