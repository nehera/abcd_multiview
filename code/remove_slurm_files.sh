# Find and remove slurm .out files
find . -name "*.out" -type f -delete

# Find and remove slurm .out files
find . -name "*.err" -type f -delete

# Find and remove jupyter .log files
find . -name "*.log" -type f -delete

# Find and remove mcmc reports
find  . -name "iter_*" -type f -delete

# BE CAREFUL WITH NEXT COMMAND
# Remove data/ outputs
rm data/*