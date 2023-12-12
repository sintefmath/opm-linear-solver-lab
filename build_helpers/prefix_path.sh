#!/bin/bash
basepath=$1

if [ $# -lt 2 ]
then
    basepath_dune=$basepath
else
    basepath_dune=$2
fi

if [ ! -d ${basepath} ]
then
    >&2 echo "$basepath does not exist (or is not a directory)"
    >&2 echo "Usage:"
    >&2 echo "    bash $0 <path_to_folder_containing_opm> <optional extra path to dune>"

    exit 1
fi

if [ ! -d ${basepath_dune} ]
then
    >&2 echo "$basepath_dune does not exist (or is not a directory)"
    >&2 echo "Usage:"
    >&2 echo "    bash $0 <path_to_folder_containing_opm> <optional extra path to dune>"

    exit 1
fi

prefix_path=""

for repo in dune-common dune-grid dune-geometry dune-istl
do
    prefix_path="$prefix_path${basepath_dune}/${repo}/build;"
done

for repo in opm-common opm-grid opm-models opm-simulators
do
    prefix_path="$prefix_path${basepath}/${repo}/build;"
done

echo $prefix_path