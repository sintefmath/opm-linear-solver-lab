import json
import os
import re


import matplotlib.pyplot as plt
import numpy as np

import plotting

def extract_solve_time(output):
    match = re.search(r"solve:\s+(\d*\.?\d*) s", output)
    return float(match.group(1)) if match else None

def extract_iterations(output):
    match = re.search(r"Iterations:\s+(\d+)", output)
    return int(match.group(1)) if match else None

def extract_error(output):
    match = re.search(r"Error:\s+([\d\.]+)", output)
    return float(match.group(1)) if match else None


def process_data(matrix_directory):
    per_iter_data = {'times': {}, 'errors': {}}

    for filename in os.listdir(matrix_directory):
        if filename.endswith('_results.json'):
            with open(os.path.join(matrix_directory, filename)) as f:
                raw_results = json.load(f)

            for method, configs_outputs in raw_results.items():
                if method not in per_iter_data['times']:
                    per_iter_data['times'][method] = {}
                    per_iter_data['errors'][method] = {}
                for config, output in configs_outputs.items():
                    solve_time = extract_solve_time(output)
                    iterations = extract_iterations(output)
                    error = extract_error(output)

                    if solve_time is not None and iterations is not None and iterations > 0:
                        time_per_iteration = (solve_time * 1e3) / iterations
                        per_iter_data['times'][method].setdefault(config, []).append(time_per_iteration)

                    if error is not None:
                        per_iter_data['errors'][method].setdefault(config, []).append(error)

    return per_iter_data


if __name__ == "__main__":
    json_files = [
        'amgcl_ilu0_options.json',
        'amgcl_ilu0_options.json -b 2',
        'amgcl_spai0_options.json',
        'amgcl_spai0_options.json -b 2',
        'amgcl_spai1_options.json',
        'amgcl_spai1_options.json -b 2'
        #'amgcl_cpr_options.json',
        #'configs/amgcl_cpr_options.json -b 2',
        #'configs/amgcl_drs_setup.json',
        #'configs/amgcl_drs_setup.json -b 2'
    ]
    CPU_commands = {
        'solver': 'CPU',
        #'cpr_drs': 'CPU (CPR DRS)'
    }

    GPU_commands = {
        'solver_cuda': 'CUDA',
        'solver_vexcl_cuda': 'VEXCL CUDA',
        #'cpr_drs_cuda': 'CUDA (CPR DRS)',
        #'cpr_drs_vexcl_cuda': 'VEXCL CUDA (CPR DRS)'
    }

    all_commands = dict(CPU_commands)
    all_commands.update(GPU_commands)

    benchmark_results_dir = 'benchmark_results/amgcl'
    subfolders = ["sleipner"]
    
    name = "CPU"
    for subfolder in subfolders:
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)  # Modify to also include error data
            plotting.plot_and_save(name, subfolder_path, subfolder, json_files, per_iter_data['times'], CPU_commands)
            plotting.plot_time_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['times'], CPU_commands)
            plotting.plot_and_save_error(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], CPU_commands)
            plotting.plot_error_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], CPU_commands)
    
    name = "GPU"
    for subfolder in subfolders:
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)  # Modify to also include error data
            plotting.plot_and_save(name, subfolder_path, subfolder, json_files, per_iter_data['times'], GPU_commands)
            plotting.plot_time_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['times'], GPU_commands)
            plotting.plot_and_save_error(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], GPU_commands)
            plotting.plot_error_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], GPU_commands)


    name = "CPU+GPU"
    for subfolder in subfolders:
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)  # Modify to also include error data
            plotting.plot_and_save(name, subfolder_path, subfolder, json_files, per_iter_data['times'], all_commands)
            plotting.plot_time_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['times'], all_commands)
            plotting.plot_and_save_error(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], all_commands)
            plotting.plot_error_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], all_commands)
