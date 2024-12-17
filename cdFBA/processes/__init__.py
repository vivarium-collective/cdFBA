from cdFBA.processes.dfba import DFBA

def register_processes(core):
    core.register_process('DFBA', DFBA)
    
    return core