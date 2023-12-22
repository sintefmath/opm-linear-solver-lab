# opm-linear-solver-lab
Experimental playground for testing out linear solvers in OPM Flow.

## Compiling
We assume you have opm and dune in your prefix path. Compiling should just be

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

Note that you probably want to specify `CMAKE_PREFIX_PATH` to point to the locations of the OPM and Dune build folder, that is

```bash
# from build. 
# IMPORTANT: Notice the quotation marks around the prefix path. 
# This is IMPORTANT when you have multiple paths
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="/path/to/opm-common/build;/path/to/opm-grid/build;MOREPATHSHERE"
```

To simplify this process, we've added the script `build_helpers/prefix_path.sh` to generate this prefix path. If you directory structure looks like this:

```bash
#DUNE
/some/path/for/dune/dune-common/build
/some/path/for/dune/dune-grid/build
/some/path/for/dune/dune-geometry/build
/some/path/for/dune/dune-istl/build

#OPM
/some/path/other/path/for/opm/opm-common/build
/some/path/other/path/for/opm/opm-grid/build
/some/path/other/path/for/opm/opm-models/build
/some/path/other/path/for/opm/opm-simulators/build
```

you can run cmake as

```bash
# from build. 
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$(bash ../build_helpers/prefix_path.sh /some/path/other/path/for/opm/ /some/path/for/dune/dune-common/)"
```

### Note on libfmt

You probably need to have libfmt installed. On Ubuntu this can be accomplished by installing

```bash
sudo apt install libfmt-dev
```

Alternatively you can manually donwload and install libfmt. Make sure to extend your `CMAKE_PREFIX_PATH` to contain the install directory of libfmt.
## Running
From the build folder, you can run eg. 

```bash
./linsolverlab \
    -x ../examples/matrices/spe1/rhs.mm \
    -y  ../examples/matrices/spe1/rhs.mm  \
    -m  ../examples/matrices/spe1/matrix.mm \
    -c ../examples/configurations/cpu/ilu0.json \
    -g ../examples/configurations/gpu/cuilu0.json
```

Example output

```JSON
{
    "CPU": {
        "runtime": "21541",
        "failed_by_exception": "false",
        "iterations": "25",
        "reduction": "9.463892568041427e-13",
        "converged": "true",
        "conv_rate": "0.33040209473507992",
        "elapsed": "0.021522673999999999",
        "condition_estimate": "-1"
    },
    "GPU": {
        "runtime": "21831",
        "failed_by_exception": "false",
        "iterations": "25",
        "reduction": "3.532348061635463e-13",
        "converged": "true",
        "conv_rate": "0.31763074967064442",
        "elapsed": "0.02177946",
        "condition_estimate": "-1"
    }
}

```


## Running only preconditioner

```bash
./linsolverlab \
    -x ../examples/matrices/spe1/rhs.mm \
    -y  ../examples/matrices/spe1/rhs.mm  \
    -m  ../examples/matrices/spe1/matrix.mm \
    -c ../examples/configurations/cpu/ilu0_onlypreconditioner.json \
    -g ../examples/configurations/gpu/cuilu0_onlypreconditioner.json
```

Example output
```JSON
{
    "CPU": {
        "runtime": "6",
        "failed_by_exception": "false",
        "iterations": "1",
        "reduction": "0",
        "converged": "true",
        "conv_rate": "1",
        "elapsed": "0",
        "condition_estimate": "-1"
    },
    "GPU": {
        "runtime": "2396",
        "failed_by_exception": "false",
        "iterations": "1",
        "reduction": "0",
        "converged": "true",
        "conv_rate": "1",
        "elapsed": "0",
        "condition_estimate": "-1"
    }
}

```


## Running benchmarks


### OPM-flow


Generating benchmark data

```python
# usage: run_opm_benchmark.py <subfolder1> [<subfolder2> ...]
python benchmarking_scripts/run_opm_benchmark.py sleipner
```

Adjusting the experiment in the file

```python
    json_files = [
        'examples/configurations/cpu/ilu0.json --block-size 2',
        'examples/configurations/cpu/dilu.json --block-size 2',
        'examples/configurations/gpu/cuilu0.json --block-size 2',
        'examples/configurations/gpu/cudilu.json --block-size 2'
    ]
    matrix_root = "./examples/matrices"  # folder where the matrix files can be found
    limit = 10  # Adjust this to change the number of matrices processed per directory
```



Plotting the results

```python
python benchmarking_scripts/process_and_plot_opm.py
```
