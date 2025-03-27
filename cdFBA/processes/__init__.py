from cdFBA.processes.dfba import dFBA, UpdateEnvironment, StaticConcentration, Injector, WaveFunction
from cdFBA.processes.dfbalauncher import EnvironmentMonitor

def register_processes(core):
    core.register_process('dFBA', dFBA)
    core.register_process('UpdateEnvironment', UpdateEnvironment)
    core.register_process('StaticConcentration', StaticConcentration)
    core.register_process('Injector', Injector)
    core.register_process('WaveFunction', WaveFunction)
    core.register_process('EnvironmentMonitor', EnvironmentMonitor)
    
    return core