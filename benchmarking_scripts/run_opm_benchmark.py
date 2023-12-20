import os
import subprocess
import json
import sys

def run_command(command):
    print(f"Running command: {command}")
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    print(f"Command completed.")
    return result.stdout

def main(binary, results_directory, json_files, matrix_root, subfolders, limit=None):
    os.makedirs(results_directory, exist_ok=True)

    # Process only specified subfolders
    for subfolder in subfolders:
        matrix_path = os.path.join(matrix_root, subfolder)
        if not os.path.isdir(matrix_path):
            print(f"Directory not found: {matrix_path}")
            continue

        matrix_files = sorted([f for f in os.listdir(matrix_path) if "matrix" in f if f.endswith('.bin')])
        rhs_files = sorted([f for f in os.listdir(matrix_path) if "rhs" in f if f.endswith('.bin')])
        
        if limit:
            matrix_files = matrix_files[:limit]
            rhs_files = rhs_files[:limit]

        for matrix_file, rhs_file in zip(matrix_files, rhs_files):
            matrix_name = os.path.splitext(matrix_file)[0]
            results = {}
            for json_file in json_files:
                config_filename = os.path.basename(json_file)
                if "cpu" in json_file:
                    label = "CPU"
                    full_cmd = f"{binary} -m '{os.path.join(matrix_path, matrix_file)}' -y '{os.path.join(matrix_path, rhs_file)}' -x '{os.path.join(matrix_path, rhs_file)}' -c {json_file}"
                elif "gpu" in json_file:
                    label = "GPU"
                    full_cmd = f"{binary} -m '{os.path.join(matrix_path, matrix_file)}' -y '{os.path.join(matrix_path, rhs_file)}' -x '{os.path.join(matrix_path, rhs_file)}' -g {json_file}"
                else:
                    print("Unrecognised config file")

                # Run the command and save the output
                output = run_command(full_cmd)
                output_dict = json.loads(output)
                if label not in results:
                    results[label] = {}
                results[label][config_filename] = output_dict[label]

                # Save the results for each matrix to a JSON file
                matrix_directory = os.path.join(results_directory, subfolder)
                os.makedirs(matrix_directory, exist_ok=True)
                with open(os.path.join(matrix_directory, f'{matrix_name}_results.json'), 'w') as f:
                    json.dump(results, f, indent=4)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_opm_benchmark.py <subfolder1> [<subfolder2> ...]")
        sys.exit(1)

    json_files = [
        'examples/configurations/cpu/dilu.json --block-size 2',
        'examples/configurations/cpu/ilu0.json --block-size 2',
        'examples/configurations/gpu/cuilu0.json --block-size 2',
        'examples/configurations/gpu/cudilu.json --block-size 2'
    ]
    matrix_root = "./examples/matrices"
    limit = 10  # Adjust this to change the number of matrices processed per directory

    binary = "./build/linsolverlab"
    results_directory = 'benchmark_results/opm'
    subfolders = sys.argv[1:]  # Get subfolder names from command-line arguments
    main(binary, results_directory, json_files, matrix_root, subfolders, limit)
