import os
import sys
import time
import numpy as np
import subprocess
from multiprocessing import Pool, current_process


BASE_DIR = os.getcwd()
SCRIPT_DIR = os.path.join(BASE_DIR, "python")
DATA_DIR = os.path.join(BASE_DIR, "data", "sim.generated")
LOG_DIR = os.path.join(BASE_DIR, "python", "pylogs")

if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)

p_mtds = ["p.somde.py", "p.spatialde2.py"]
r_mtds = ["r.somde.py", "r.spatialde2.py"]
i_mtds = ["i.somde.py", "i.spatialde2.py"] 

si = np.ones(1)

N_CORES = 1

def run_python_method(args):
    """run python method in subprocess with logging."""
    method, siv, sim_type = args
    process_id = current_process().pid
    
    print(f">>> [Worker {process_id}] Processing {sim_type.upper()} SI: {siv:.2f} with {method}")
    
    cmd = f"taskset -c 0-11 python {os.path.join(SCRIPT_DIR, method)} {siv}"
    
    method_name = method.replace('.py', '').replace('.', '_')
    log_file = os.path.join(LOG_DIR, f"{method_name}.{sim_type}.SI{siv:.2f}.log")
    
    cmd_with_log = f"{cmd} > {log_file} 2>&1"
    
    try:
        exit_code = subprocess.run(cmd_with_log, shell=True, check=False).returncode
        if exit_code != 0:
            error_msg = f"!!! [Worker {process_id}] Error in method {method} at SI {siv:.2f}"
            print(error_msg)
            return f"ERROR: {method} SI {siv:.2f}"
        else:
            return f"SUCCESS: {method} SI {siv:.2f}"
    except Exception as e:
        error_msg = f"!!! [Worker {process_id}] Exception in {method} at SI {siv:.2f}: {str(e)}"
        print(error_msg)
        return f"EXCEPTION: {method} SI {siv:.2f}"

if __name__ == "__main__":
    t_start = time.time()
    
    print("Preparing parallel simulations...")

    print("\n=== Starting progression benchmark simulations... ===")
    with Pool(N_CORES) as pool:
        p_tasks = [(method, siv, "p") for siv in si for method in p_mtds]
        p_results = pool.map(run_python_method, p_tasks)
    
    print("\n Progression Done.")
    
    print("\n=== Starting regionation benchmark simulations... ===")
    with Pool(N_CORES) as pool:
        r_tasks = [(method, siv, "r") for siv in si for method in r_mtds]
        r_results = pool.map(run_python_method, r_tasks)
    
    print("\n Regionation Done.")
    
    # Islet 시뮬레이션  
    print("\n=== Starting islet benchmark simulations... ===")
    with Pool(N_CORES) as pool:
        i_tasks = [(method, siv, "i") for siv in si for method in i_mtds]
        i_results = pool.map(run_python_method, i_tasks)
    
    print("\n Islet Done.")
    print("\n=== All simulations completed! ===")
    
    elapsed_hours = (time.time() - t_start) / 3600
    print(f"\n Elapsed time: {elapsed_hours:.2f} hours")