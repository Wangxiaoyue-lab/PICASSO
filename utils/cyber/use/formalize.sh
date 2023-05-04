#!/bin/bash

task=$1
dir=$2

mkdir ${task}

if [ ! -d ${task}/script ]; then
    mkdir ${task}/script
fi

if [ ! -d ${task}/input ]; then
    mkdir ${task}/input
fi

if [ ! -d ${task}/output ]; then
    mkdir ${task}/output
fi

if [ ! -d ${task}/report ]; then
    mkdir ${task}/report
fi

if [ ! -d ${task}/output/picture ]; then
    mkdir ${task}/output/picture
fi

if [ ! -d ${task}/output/table ]; then
    mkdir ${task}/output/table
fi

if [ ! -f ${task}/input/description.yaml ]; then
    touch ${task}/input/description.yaml
fi

echo "Input1:\n  Name:\n  Date:\n  Source:\nInput2:\n  Name:\n  Date:\n  Source:" > ${task}/input/description.yaml

if [ ! -f ${task}/report/report.yaml ]; then
    touch ${task}/report/report.yaml
fi
echo "Experiment_name: _name_\n:Experiment_situation:\n  Experiment_location:_location_\n  Experiment_operater: _operater_\n  Experiment_date: _date_\n  Experiment_platform: _platform_\nExperiment_aim: \nexperiment_methods_and_material:\n  Experiment_data:_data_\n  Experiment_software:\n    - _software_\n  Supplement: _Supplement_\nExperiment_process: _process_\nExperiment_result: \n  Experiment_output: _output_\n  Experiment_error: _error_\nExperiment_discussion: _discussion_\nNext_plan: _plan_" >  ${task}/report/report.yaml