#!/bin/bash

# Define the number of tasks and memory limit
num_tasks_parallel=6  # Adjust based on your CPU cores
num_tasks_overall=9 # Adjust based on the ntasks to execute
memory_limit=16g  # Adjust based on your local machine memory

# Change directory to the location of your R script
cd ~/Documents/GitHub/abcd_multiview

# Define the list of task IDs
task_ids=$(seq 1 $num_tasks_overall)

# Run the R script in parallel with different task IDs
parallel --jobs $num_tasks_parallel Rscript simulation_study.R {} ::: $task_ids
