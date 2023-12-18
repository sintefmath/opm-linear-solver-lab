import json
import os
import re


import matplotlib.pyplot as plt
import numpy as np


def extract_solve_time(output):
    match = re.search(r"solve:\s+(\d*\.?\d*) s", output)
    return float(match.group(1)) if match else None

def extract_iterations(output):
    match = re.search(r"Iterations:\s+(\d+)", output)
    return int(match.group(1)) if match else None

def extract_error(output):
    match = re.search(r"Error:\s+([\d\.]+)", output)
    return float(match.group(1)) if match else None

def plot_and_save(name, matrix, json_files, per_iter_times, commands):
    bar_width = 0.2
    num_configs = len(json_files)
    positions = np.arange(len(commands))

    # Define a list of colors, one for each configuration
    colors = ['blue', 'green', 'red', 'orange', 'purple',  'cyan', 'pink', 'yellow', 'grey', 'brown']
    if len(json_files) > len(colors):
        print("Warning: Not enough colors defined for the number of configurations.")

    plt.figure(figsize=(14, 6))

    # Used to track if the label has been added to the legend
    added_labels = set()

    for i, config_name in enumerate(json_files):
        config_filename = os.path.basename(config_name)
        color = colors[i % len(colors)]  # Use modulo to cycle through colors if necessary
        for j, method in enumerate(commands.values()):
            times = per_iter_times[method].get(config_filename, [])
            if times:
                avg_time = np.mean(times)
                std_dev = np.std(times)
                label = config_filename if config_filename not in added_labels else None
                plt.bar(positions[j] + i * bar_width, avg_time, yerr=std_dev, width=bar_width, label=label, color=color)
                added_labels.add(config_filename)  # Mark this label as added

    plt.xlabel('Method')
    plt.ylabel('Average Time per Iteration (ms) with Std Dev')
    plt.title(f'{name}: Average Time Per Iteration for {matrix} Matrix with Std Dev')
    plt.xticks(positions + bar_width * (num_configs / 2), list(commands.values()))

    plt.legend()
    matrix_directory = os.path.join('benchmark_results', matrix)
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_time_per_iteration.png"))
    plt.close()

def plot_and_save_error(name, matrix, json_files, per_iter_errors, commands):
    bar_width = 0.2
    num_configs = len(json_files)
    positions = np.arange(len(commands))

    # Define a list of colors, one for each configuration
    colors = ['blue', 'green', 'red', 'orange', 'purple', 'cyan', 'pink', 'yellow', 'grey', 'brown']
    if len(json_files) > len(colors):
        print("Warning: Not enough colors defined for the number of configurations.")

    plt.figure(figsize=(14, 6))

    # Used to track if the label has been added to the legend
    added_labels = set()

    for i, config_name in enumerate(json_files):
        config_filename = os.path.basename(config_name)
        color = colors[i % len(colors)]  # Use modulo to cycle through colors if necessary
        for j, method in enumerate(commands.values()):
            errors = per_iter_errors[method].get(config_filename, [])
            if errors:
                avg_error = np.mean(errors)
                std_dev = np.std(errors)
                label = config_filename if config_filename not in added_labels else None
                plt.bar(positions[j] + i * bar_width, avg_error, yerr=std_dev, width=bar_width, label=label, color=color)
                added_labels.add(config_filename)  # Mark this label as added

    plt.xlabel('Method')
    plt.ylabel('Average Error with Std Dev')
    plt.title(f'{name}: Average Error for {matrix} Matrix with Std Dev')
    plt.xticks(positions + bar_width * (num_configs / 2), list(commands.values()))
    plt.yscale('log')  # Set y-axis to logarithmic scale
    plt.legend()
    matrix_directory = os.path.join('benchmark_results', matrix)
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_error.png"))
    plt.close()



def plot_time_histogram(name, matrix, per_iter_times, commands):
    plt.figure(figsize=(14, 12))  # Adjust the figure size as needed

    for i, method in enumerate(commands.values(), start=1):
        plt.subplot(2, 3, i)  # Arrange subplots in a 2x3 grid
        times = [time for config_times in per_iter_times.get(method, {}).values() for time in config_times]
        if times:
            plt.hist(times, 20, alpha=0.5)
            plt.xlabel('Time per Iteration (ms)')
            plt.ylabel('Frequency')
            plt.title(f'{method}')

    matrix_directory = os.path.join('benchmark_results', matrix)
    plt.suptitle(f'Histogram of Time Per Iteration for Matrices in {os.path.basename(matrix_directory)}')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_time_per_iteration_histogram.png"))
    plt.close()



def plot_error_histogram(name, matrix, per_iter_errors, commands):
    plt.figure(figsize=(14, 12))  # Adjust the figure size as needed

    for i, method in enumerate(commands.values(), start=1):
        plt.subplot(2, 3, i)  # Arrange subplots in a 2x3 grid
        errors = [error for config_errors in per_iter_errors.get(method, {}).values() for error in config_errors]
        if errors:
            plt.hist(errors, 20, alpha=0.5)
            plt.xlabel('Error')
            plt.ylabel('Frequency')
            plt.title(f'{method}')
            #plt.xscale('log')  # Set y-axis to logarithmic scale

    matrix_directory = os.path.join('benchmark_results', matrix)
    plt.suptitle(f'{name}: Histogram of Error for Matrices in {os.path.basename(matrix_directory)}')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_error_histogram.png"))
    plt.close()

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

    benchmark_results_dir = 'benchmark_results'
    
    name = "CPU"
    for subfolder in os.listdir(benchmark_results_dir):
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)  # Modify to also include error data
            plot_and_save(name, subfolder, json_files, per_iter_data['times'], CPU_commands)
            plot_time_histogram(name, subfolder, per_iter_data['times'], CPU_commands)
            plot_and_save_error(name, subfolder, json_files, per_iter_data['errors'], CPU_commands)
            plot_error_histogram(name, subfolder, per_iter_data['errors'], CPU_commands)

    name = "GPU"
    for subfolder in os.listdir(benchmark_results_dir):
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)  # Modify to also include error data
            plot_and_save(name, subfolder, json_files, per_iter_data['times'], GPU_commands)
            plot_time_histogram(name, subfolder, per_iter_data['times'], GPU_commands)
            plot_and_save_error(name, subfolder, json_files, per_iter_data['errors'], GPU_commands)
            plot_error_histogram(name, subfolder, per_iter_data['errors'], GPU_commands)


    name = "CPU+GPU"
    for subfolder in os.listdir(benchmark_results_dir):
        subfolder_path = os.path.join(benchmark_results_dir, subfolder)
        if os.path.isdir(subfolder_path):
            per_iter_data = process_data(subfolder_path)  # Modify to also include error data
            plot_and_save(name, subfolder, json_files, per_iter_data['times'], all_commands)
            plot_time_histogram(name, subfolder, per_iter_data['times'], all_commands)
            plot_and_save_error(name, subfolder, json_files, per_iter_data['errors'], all_commands)
            plot_error_histogram(name, subfolder, per_iter_data['errors'], all_commands)