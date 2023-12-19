import json
import os


import matplotlib.pyplot as plt
import numpy as np



def plot_and_save(name, matrix_directory, matrix, json_files, per_iter_times, commands):    

    bar_width = 0.2
    positions = np.arange(len(commands))
    plt.figure(figsize=(14, 6))

    labels = []
    ticks = []

    for i, method in enumerate(commands.values()):
            config_ticks = []
            for j, config_name in enumerate(json_files):
                config_filename = os.path.basename(config_name)
                times = per_iter_times[method].get(config_filename, [])
                if times:
                    avg_time = np.mean(times)
                    std_dev = np.std(times)
                    plt.bar(positions[i] + j * bar_width, avg_time, yerr=std_dev, width=bar_width, label=config_filename)#, color=color)
                    config_ticks.append(positions[i] + j * bar_width)
            ticks.append(np.mean(config_ticks))
            labels.append(method)

    plt.ylabel('Average Time per Iteration (ms) with Std Dev')
    plt.title(f'{name}: Average Time Per Iteration for {matrix} Matrix with Std Dev')
    
    if len(commands) > 1:
        plt.xticks(ticks, labels)
        plt.xlabel('Method')
    else:
        plt.xticks([])
    
    plt.legend()
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_time_per_iteration.png"))
    plt.close()

def plot_and_save_error(name, matrix_directory, matrix, json_files, per_iter_errors, commands):
    bar_width = 0.2
    num_configs = len(json_files)
    positions = np.arange(len(commands))

    plt.figure(figsize=(14, 6))

    labels = []
    ticks = []
    for i, method in enumerate(commands.values()):
            config_ticks = []
            for j, config_name in enumerate(json_files):
                config_filename = os.path.basename(config_name)
                errors = per_iter_errors[method].get(config_filename, [])
                if errors:
                    avg_error = np.mean(errors)
                    std_dev = np.std(errors)
                    plt.bar(positions[i] + j * bar_width, avg_error, yerr=std_dev, width=bar_width, label=config_filename)#, color=color)
                    config_ticks.append(positions[i] + j * bar_width)
            ticks.append(np.mean(config_ticks))
            labels.append(method)

    plt.xlabel('Method')
    plt.ylabel('Average Error with Std Dev')
    plt.title(f'{name}: Average Error for {matrix} Matrix with Std Dev')
    plt.xticks(positions + bar_width * (num_configs / 2), list(commands.values()))
    plt.yscale('log')  # Set y-axis to logarithmic scale
    plt.legend()
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_error.png"))
    plt.close()

def plot_time_histogram(name, matrix_directory, matrix, json_files, per_iter_times, commands):
    num_methods = len(commands)
    num_configs = len(json_files)
    total_plots = num_methods * num_configs

    # Calculate the number of rows and columns for the subplots
    # assuming a maximum of 3 columns
    num_cols = min(3, total_plots)
    num_rows = (total_plots + num_cols - 1) // num_cols  # Ceiling division

    plt.figure(figsize=(14, 6 * num_rows))  # Adjust the figure size as needed

    subplot_index = 1
    for method in commands.values():
        for config_name in json_files:
            config_filename = os.path.basename(config_name)
            
            times = per_iter_times[method].get(config_filename, [])
            if times:
                plt.subplot(num_rows, num_cols, subplot_index)
                plt.hist(times, 20, alpha=0.5)
                plt.xlabel('Time per Iteration (ms)')
                plt.ylabel('Frequency')
                plt.title(f'{method} - {config_filename}')

                subplot_index += 1

    plt.suptitle(f'Histogram of Time Per Iteration for Matrices in {matrix}')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_time_per_iteration_histogram.png"))
    plt.close()


def plot_error_histogram(name, matrix_directory, matrix, json_files, per_iter_errors, commands):
    num_methods = len(commands)
    num_configs = len(json_files)
    total_plots = num_methods * num_configs

    # Calculate the number of rows and columns for the subplots
    # assuming a maximum of 3 columns
    num_cols = min(3, total_plots)
    num_rows = (total_plots + num_cols - 1) // num_cols  # Ceiling division

    plt.figure(figsize=(14, 6 * num_rows))  # Adjust the figure size as needed

    subplot_index = 1
    for method in commands.values():
        for config_name in json_files:
            config_filename = os.path.basename(config_name)
            
            times = per_iter_errors[method].get(config_filename, [])
            if times:
                plt.subplot(num_rows, num_cols, subplot_index)
                plt.hist(times, 20, alpha=0.5)
                plt.xlabel('Error')
                plt.ylabel('Frequency')
                plt.title(f'{method} - {config_filename}')

                subplot_index += 1

    plt.suptitle(f'Histogram of Error for Matrices in {matrix}')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_error_histogram.png"))
    plt.close()

def process_data(matrix_directory):
    per_iter_data = {'times': {}, 'errors': {}}

    for filename in os.listdir(matrix_directory):
        if filename.endswith('_results.json'):
            with open(os.path.join(matrix_directory, filename)) as f:
                raw_results = json.load(f)

            for method in raw_results:
                if method not in per_iter_data['times']:
                    per_iter_data['times'][method] = {}
                    per_iter_data['errors'][method] = {}

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

                    if not np.isnan(reduction):
                        per_iter_data['errors'][method].setdefault(config, []).append(reduction)

    return per_iter_data