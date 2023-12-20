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
                    plt.bar(positions[i] + j * bar_width, avg_time, yerr=std_dev, width=bar_width, label=config_filename)
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


def plot_and_save_solve(name, matrix_directory, matrix, json_files, solve_times, commands):

    bar_width = 0.2
    positions = np.arange(len(commands))
    plt.figure(figsize=(14, 6))

    labels = []
    ticks = []

    for i, method in enumerate(commands.values()):
            config_ticks = []
            for j, config_name in enumerate(json_files):
                config_filename = os.path.basename(config_name)
                solve_time = solve_times[method].get(config_filename, [])
                if solve_time:
                    solve_time = np.mean(solve_time)
                    std_dev = np.std(solve_time)
                    plt.bar(positions[i] + j * bar_width, solve_time, yerr=std_dev, width=bar_width, label=config_filename)
                    config_ticks.append(positions[i] + j * bar_width)
            ticks.append(np.mean(config_ticks))
            labels.append(method)

    plt.ylabel('Average Solve Time (s) with Std Dev')
    plt.title(f'{name}: Average Solve Time for {matrix} Matrix with Std Dev')

    if len(commands) > 1:
        plt.xticks(ticks, labels)
        plt.xlabel('Method')
    else:
        plt.xticks([])

    plt.legend()
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_solve_time.png"))
    plt.close()


def plot_and_save_iter(name, matrix_directory, matrix, json_files, iterations, commands):    

    bar_width = 0.2
    positions = np.arange(len(commands))
    plt.figure(figsize=(14, 6))

    labels = []
    ticks = []

    for i, method in enumerate(commands.values()):
            config_ticks = []
            for j, config_name in enumerate(json_files):
                config_filename = os.path.basename(config_name)
                iter = iterations[method].get(config_filename, [])
                if iter:
                    iter_time = np.mean(iter)
                    std_dev = np.std(iter)
                    plt.bar(positions[i] + j * bar_width, iter_time, yerr=std_dev, width=bar_width, label=config_filename)
                    config_ticks.append(positions[i] + j * bar_width)
            ticks.append(np.mean(config_ticks))
            labels.append(method)

    plt.ylabel('Average Number of Iterations with Std Dev')
    plt.title(f'{name}: Average Number of iterations for {matrix} Matrix with Std Dev')

    if len(commands) > 1:
        plt.xticks(ticks, labels)
        plt.xlabel('Method')
    else:
        plt.xticks([])
    plt.legend()
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_iterations.png"))
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
                    plt.bar(positions[i] + j * bar_width, avg_error, yerr=std_dev, width=bar_width, label=config_filename)
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
