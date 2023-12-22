import os

import matplotlib.pyplot as plt
import numpy as np


def plot_and_save_time_per_iter(name, matrix_directory, matrix, json_files, per_iter_times, commands):
    
    bar_width = 0.2
    positions = np.arange(len(commands))
    fig, ax = plt.subplots(figsize=(14, 6))

    # Create a color map for each config_filename
    color_map = plt.cm.get_cmap('tab10')
    colors = {config_name: color_map(i) for i, config_name in enumerate(json_files)}

    method_labels = []
    method_ticks = []
    bars_added = 0
    for i, method in enumerate(commands.values()):
        config_ticks = []
        if i > 0:
            bars_added = 1  # add spacing between methods
        for j, config_name in enumerate(json_files):
            config_filename = os.path.basename(config_name)
            times = per_iter_times[method].get(config_filename, [])
            if times:
                avg_time = np.mean(times)
                std_dev = np.std(times)
                bar_color = colors[config_name]
                handles, labels = ax.get_legend_handles_labels()
                if not config_filename in labels:
                    label = config_filename
                else:
                    label = ""
                bar = ax.bar(positions[i] + bars_added * bar_width, avg_time, yerr=std_dev, width=bar_width, color=bar_color, label=label, capsize=5)
                config_ticks.append(positions[i] + bars_added * bar_width)
                bars_added += 1
        method_ticks.append(np.mean(config_ticks))
        method_labels.append(method)
    
    ax.set_ylabel('Average Time per Iteration (ms) with Std Dev', fontsize=12)
    ax.set_title(f'{name}: Average Time Per Iteration for {matrix} Matrix with Std Dev', fontsize=14)
    if len(commands) > 1:
        ax.set_xticks(method_ticks)
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.set_xlabel('Method', fontsize=12)
    else:
        ax.set_xticks([])
    ax.legend()
    ax.grid(True, which='major', linestyle='--', linewidth='0.5', color='grey')
    plt.tight_layout() 
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_time_per_iteration.png"), dpi=300)
    plt.close()


def plot_and_save_solve(name, matrix_directory, matrix, json_files, solve_times, commands):

    bar_width = 0.2
    positions = np.arange(len(commands))
    fig, ax = plt.subplots(figsize=(14, 6))

    # Create a color map for each config_filename
    color_map = plt.cm.get_cmap('tab10')
    colors = {config_name: color_map(i) for i, config_name in enumerate(json_files)}

    method_labels = []
    method_ticks = []
    bars_added = 0
    for i, method in enumerate(commands.values()):
        config_ticks = []
        if i > 0:
            bars_added = 1  # add spacing between methods
        for j, config_name in enumerate(json_files):
            config_filename = os.path.basename(config_name)
            solve_time = solve_times[method].get(config_filename, [])
            if solve_time:
                avg_solve_time = np.mean(solve_time)
                std_dev = np.std(solve_time)
                bar_color = colors[config_name]
                handles, labels = ax.get_legend_handles_labels()
                if not config_filename in labels:
                    label = config_filename
                else:
                    label = ""
                ax.bar(positions[i] + bars_added * bar_width, avg_solve_time, yerr=std_dev, width=bar_width, color=bar_color, label=label, capsize=5)
                config_ticks.append(positions[i] + bars_added * bar_width)
                bars_added += 1
        method_ticks.append(np.mean(config_ticks))
        method_labels.append(method)

    ax.set_ylabel('Average Solve Time (s) with Std Dev', fontsize=12)
    ax.set_title(f'{name}: Average Solve Time for {matrix} Matrix with Std Dev', fontsize=14)

    if len(commands) > 1:
        ax.set_xticks(method_ticks)
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.set_xlabel('Method', fontsize=12)
    else:
        ax.set_xticks([])

    ax.legend()
    ax.grid(True, which='major', linestyle='--', linewidth='0.5', color='grey')
    plt.tight_layout()
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_solve_time.png"), dpi=300)
    plt.close()

def plot_and_save_iter(name, matrix_directory, matrix, json_files, iterations, commands):

    bar_width = 0.2
    positions = np.arange(len(commands))
    fig, ax = plt.subplots(figsize=(14, 6))

    # Create a color map for each config_filename
    color_map = plt.cm.get_cmap('tab10')
    colors = {config_name: color_map(i) for i, config_name in enumerate(json_files)}

    method_labels = []
    method_ticks = []
    bars_added = 0
    for i, method in enumerate(commands.values()):
        config_ticks = []
        if i > 0:
            bars_added = 1  # add spacing between methods
        for j, config_name in enumerate(json_files):
            config_filename = os.path.basename(config_name)
            iteration = iterations[method].get(config_filename, [])
            if iteration:
                iter_time = np.mean(iteration)
                std_dev = np.std(iteration)
                bar_color = colors[config_name]
                handles, labels = ax.get_legend_handles_labels()
                if not config_filename in labels:
                    label = config_filename
                else:
                    label = ""
                ax.bar(positions[i] + bars_added * bar_width, iter_time, yerr=std_dev, width=bar_width, color=bar_color, label=label, capsize=5)
                config_ticks.append(positions[i] + bars_added * bar_width)
                bars_added += 1
        method_ticks.append(np.mean(config_ticks))
        method_labels.append(method)

    ax.set_ylabel('Average Number of Iterations with Std Dev', fontsize=12)
    ax.set_title(f'{name}: Average Number of Iterations for {matrix} Matrix with Std Dev', fontsize=14)

    if len(commands) > 1:
        ax.set_xticks(method_ticks)
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.set_xlabel('Method', fontsize=12)
    else:
        ax.set_xticks([])
    ax.legend()
    ax.grid(True, which='major', linestyle='--', linewidth='0.5', color='grey')
    plt.tight_layout()
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_iterations.png"), dpi=300)
    plt.close()

def plot_and_save_error(name, matrix_directory, matrix, json_files, per_iter_errors, commands):
    bar_width = 0.2
    positions = np.arange(len(commands))
    fig, ax = plt.subplots(figsize=(14, 6))

    # Create a color map for each config_filename
    color_map = plt.cm.get_cmap('tab10')
    colors = {config_name: color_map(i) for i, config_name in enumerate(json_files)}

    method_labels = []
    method_ticks = []
    bars_added = 0
    for i, method in enumerate(commands.values()):
        config_ticks = []
        if i > 0:
            bars_added = 1  # add spacing between methods
        for j, config_name in enumerate(json_files):
            config_filename = os.path.basename(config_name)
            errors = per_iter_errors[method].get(config_filename, [])
            if errors:
                avg_error = np.mean(errors)
                std_dev = np.std(errors)
                bar_color = colors[config_name]
                handles, labels = ax.get_legend_handles_labels()
                if not config_filename in labels:
                    label = config_filename
                else:
                    label = ""
                ax.bar(positions[i] + bars_added * bar_width, avg_error, yerr=std_dev, width=bar_width, color=bar_color, label=label, capsize=5)
                config_ticks.append(positions[i] + bars_added * bar_width)
                bars_added += 1
        method_ticks.append(np.mean(config_ticks))
        method_labels.append(method)

    if len(commands) > 1:
        ax.set_xticks(method_ticks)
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.set_xlabel('Method', fontsize=12)
    else:
        ax.set_xticks([])

    ax.set_ylabel('Average Error with Std Dev', fontsize=12)
    ax.set_title(f'{name}: Average Error for {matrix} Matrix with Std Dev', fontsize=14)
    ax.set_yscale('log')  # Set y-axis to logarithmic scale
    ax.legend()
    ax.grid(True, which='major', linestyle='--', linewidth='0.5', color='grey')
    plt.tight_layout()
    plt.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_average_error.png"), dpi=300)
    plt.close()

def plot_time_histogram(name, matrix_directory, matrix, json_files, per_iter_times, commands):
    total_plots = sum(bool(per_iter_times[method].get(os.path.basename(config_name), [])) for method in commands.values() for config_name in json_files)
    
    # Calculate the number of rows and columns for the subplots
    num_cols = min(3, total_plots)
    num_rows = (total_plots + num_cols - 1) // num_cols  # Ceiling division
    
    # Create subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 4 * num_rows), squeeze=False)  # Prevents axes from being squeezed into a 1D array
    fig.suptitle(f'Histogram of Time Per Iteration for Matrices in {matrix}', fontsize=16)

    # Flatten the axes array for easy iteration, regardless of its dimensions
    axes = axes.flatten()
    
    subplot_index = 0
    for method in commands.values():
        for config_name in json_files:
            config_filename = os.path.basename(config_name)
            times = per_iter_times[method].get(config_filename, [])
            if times:
                ax = axes[subplot_index]
                ax.hist(times, bins=20, alpha=0.7, color='cornflowerblue', edgecolor='black')
                ax.set_xlabel('Time per Iteration (ms)', fontsize=12)
                ax.set_ylabel('Frequency', fontsize=12)
                ax.set_title(f'{method} - {config_filename}', fontsize=14)
                ax.grid(True, which='major', linestyle='--', linewidth='0.5', color='grey')
                ax.set_xlim([min(times), max(times)])  # Optional: Set to standardize across all plots
                subplot_index += 1

    # Turn off any empty subplots
    for idx in range(subplot_index, num_rows * num_cols):
        axes[idx].axis('off')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust the layout
    fig.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_time_per_iteration_histogram.png"), dpi=300)
    plt.close(fig)

def plot_error_histogram(name, matrix_directory, matrix, json_files, per_iter_errors, commands):
    total_plots = sum(bool(per_iter_errors[method].get(os.path.basename(config_name), [])) for method in commands.values() for config_name in json_files)

    # Calculate the number of rows and columns for the subplots
    num_cols = min(3, total_plots)
    num_rows = (total_plots + num_cols - 1) // num_cols  # Ceiling division

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 4 * num_rows), squeeze=False)
    fig.suptitle(f'Histogram of Error for Matrices in {matrix}', fontsize=16)

    # Flatten the axes array for easy iteration
    axes = axes.flatten()

    subplot_index = 0
    for method in commands.values():
        for config_name in json_files:
            config_filename = os.path.basename(config_name)
            errors = per_iter_errors[method].get(config_filename, [])
            if errors:
                ax = axes[subplot_index]

                # Define custom bins with more focus on the smaller values
                max_error = max(errors)
                bins = np.logspace(np.log10(min(errors)), np.log10(max_error), 50)
                
                # Plot histogram with logarithmic scale
                ax.hist(errors, bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
                ax.set_xscale('log')  # Set x-axis to logarithmic scale
                ax.set_xlabel('Error (log scale)', fontsize=12)
                ax.set_ylabel('Frequency', fontsize=12)
                ax.set_title(f'{method} - {config_filename}', fontsize=14)
                ax.grid(True, which='major', linestyle='--', linewidth='0.5', color='grey')
                
                subplot_index += 1

    # Turn off any empty subplots
    for idx in range(subplot_index, num_rows * num_cols):
        axes[idx].axis('off')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(os.path.join(matrix_directory, f"{name}_{matrix}_error_histogram.png"), dpi=300)
    plt.close(fig)
