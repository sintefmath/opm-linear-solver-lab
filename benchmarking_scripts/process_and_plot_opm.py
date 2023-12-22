import json
import os

import numpy as np

import plotting


def process_data(matrix_directory):
    per_iter_data = {'times': {}, 'errors': {}, 'solve_time': {}, 'iterations': {}}

    for filename in os.listdir(matrix_directory):
        if filename.endswith('_results.json'):
            with open(os.path.join(matrix_directory, filename)) as f:
                raw_results = json.load(f)

            for method in raw_results:
                if method not in per_iter_data['times']:
                    per_iter_data['times'][method] = {}
                    per_iter_data['errors'][method] = {}
                    per_iter_data['solve_time'][method] = {}
                    per_iter_data['iterations'][method] = {}

                for config in raw_results[method]:
                    output = raw_results[method][config]
                    # Convert string to boolean
                    failed_by_exception = output["failed_by_exception"] == "true"

                    # Convert string to float or int
                    iterations = int(output["iterations"])
                    elapsed = float(output["elapsed"])
                    reduction = output["reduction"]
                    reduction = float(reduction) if reduction not in ["-nan", "nan"] else np.nan

                    if not failed_by_exception and iterations > 0:
                        time_per_iteration = (elapsed * 1e3) / iterations
                        per_iter_data['times'][method].setdefault(config, []).append(time_per_iteration)
                        per_iter_data['solve_time'][method].setdefault(config, []).append(elapsed)
                        per_iter_data['iterations'][method].setdefault(config, []).append(iterations)

                    if not np.isnan(reduction):
                        per_iter_data['errors'][method].setdefault(config, []).append(reduction)

    return per_iter_data


if __name__ == "__main__":
    json_files = [
        'examples/configurations/cpu/ilu0.json --block-size 2',
        'examples/configurations/cpu/dilu.json --block-size 2',
        'examples/configurations/gpu/cuilu0.json --block-size 2',
        'examples/configurations/gpu/cudilu.json --block-size 2'
    ]

    CPU_commands = {
        'CPU': 'CPU'
    }

    GPU_commands = {
        'GPU': "GPU"
    }

    all_commands = dict(CPU_commands)
    all_commands.update(GPU_commands)

    benchmark_results_dir = 'benchmark_results/opm'
    subfolders = ["sleipner", "refined_sleipner"]

    name = "CPU"
    for subfolder in subfolders:
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)
            plotting.plot_and_save_time_per_iter(name, subfolder_path, subfolder, json_files, per_iter_data['times'], CPU_commands)
            plotting.plot_and_save_solve(name, subfolder_path, subfolder, json_files, per_iter_data['solve_time'], CPU_commands)
            plotting.plot_and_save_iter(name, subfolder_path, subfolder, json_files, per_iter_data['iterations'], CPU_commands)
            plotting.plot_time_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['times'], CPU_commands)
            plotting.plot_and_save_error(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], CPU_commands)
            plotting.plot_error_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], CPU_commands)

    name = "GPU"
    for subfolder in subfolders:
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)
            plotting.plot_and_save_time_per_iter(name, subfolder_path, subfolder, json_files, per_iter_data['times'], GPU_commands)
            plotting.plot_and_save_solve(name, subfolder_path, subfolder, json_files, per_iter_data['solve_time'], GPU_commands)
            plotting.plot_and_save_iter(name, subfolder_path, subfolder, json_files, per_iter_data['iterations'], GPU_commands)
            plotting.plot_time_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['times'], GPU_commands)
            plotting.plot_and_save_error(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], GPU_commands)
            plotting.plot_error_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], GPU_commands)


    name = "CPU+GPU"
    for subfolder in subfolders:
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)
            plotting.plot_and_save_time_per_iter(name, subfolder_path, subfolder, json_files, per_iter_data['times'], all_commands)
            plotting.plot_and_save_solve(name, subfolder_path, subfolder, json_files, per_iter_data['solve_time'], all_commands)
            plotting.plot_and_save_iter(name, subfolder_path, subfolder, json_files, per_iter_data['iterations'], all_commands)
            plotting.plot_time_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['times'], all_commands)
            plotting.plot_and_save_error(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], all_commands)
            plotting.plot_error_histogram(name, subfolder_path, subfolder, json_files, per_iter_data['errors'], all_commands)
