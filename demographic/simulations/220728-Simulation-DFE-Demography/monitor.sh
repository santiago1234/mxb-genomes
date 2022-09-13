echo $(ls results/simulations/sim-{1..350}-pop.bin 2>&1 >/dev/null  |wc -l) missing
echo $(find results/simulations/ |wc -l) completed
