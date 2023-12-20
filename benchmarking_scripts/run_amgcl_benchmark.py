import os
import subprocess
import json
import sys

def run_command(command):
    print(f"Running command: {command}")
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    print(f"Command completed.")
    return result.stdout

def main(binary, results_directory, json_files, matrix_root, commands, subfolders, limit=None):
    os.makedirs(results_directory, exist_ok=True)

    # Process only specified subfolders
    for subfolder in subfolders:
        matrix_path = os.path.join(matrix_root, subfolder)
        #matrix_path = os.path.join(matrix_path, "reports")
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
                for cmd_key, label in commands.items():
                    if 'cpr' in cmd_key and 'cpr' not in config_filename:
                        continue
                    if 'cpr' not in cmd_key and 'cpr' in config_filename:
                        continue

                    full_cmd = f"{binary}/{cmd_key} -B -A '{os.path.join(matrix_path, matrix_file)}' -f '{os.path.join(matrix_path, rhs_file)}' -P {json_file}"
                    
                    # Run the command and save the output
                    output = run_command(full_cmd)
                    if label not in results:
                        results[label] = {}
                    results[label][config_filename] = output

                    # Save the results for each matrix to a JSON file
                    matrix_directory = os.path.join(results_directory, subfolder)
                    os.makedirs(matrix_directory, exist_ok=True)
                    with open(os.path.join(matrix_directory, f'{matrix_name}_results.json'), 'w') as f:
                        json.dump(results, f, indent=4)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <subfolder1> [<subfolder2> ...]")
        sys.exit(1)

    json_files = [
        'examples/configurations/amgcl/amgcl_ilu0_options.json',
        'examples/configurations/amgcl/amgcl_ilu0_options.json -b 2',
        'examples/configurations/amgcl/amgcl_spai0_options.json',
        'examples/configurations/amgcl/amgcl_spai0_options.json -b 2',
        'examples/configurations/amgcl/amgcl_spai1_options.json',
        'examples/configurations/amgcl/amgcl_spai1_options.json -b 2',
        'examples/configurations/amgcl/amgcl_cpr_options.json',
        'examples/configurations/amgcl/amgcl_cpr_options.json -b 2',
        'examples/configurations/amgcl/amgcl_cpr_drs_setup.json',
        'examples/configurations/amgcl/amgcl_cpr_drs_setup.json -b 2'
    ]
    matrix_root = "examples/matrices"
    commands = {
        'solver': 'CPU',
        'solver_cuda': 'CUDA',
        #'solver_vexcl_cuda': 'VEXCL CUDA',
        #'cpr_drs': 'CPU (CPR DRS)',
        #'cpr_drs_cuda': 'CUDA (CPR DRS)',
        #'cpr_drs_vexcl_cuda': 'VEXCL CUDA (CPR DRS)'
    }
    limit = 10  # Adjust this to change the number of matrices processed per directory

    binary = "/home/jakob/code/amgcl/build/examples"
    results_directory = 'benchmark_results/amgcl'
    subfolders = sys.argv[1:]  # Get subfolder names from command-line arguments
    main(binary, results_directory, json_files, matrix_root, commands, subfolders, limit)
