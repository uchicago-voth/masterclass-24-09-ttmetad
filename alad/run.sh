module load gromacs/2021.5
gmx_mpi mdrun -v -deffnm step5_production -ntomp 4 -plumed plumed.dat

# plumed sum_hills --hills HILLS --mintozero --stride 100 --outfile fes/fes_ttmetad
# plumed sum_hills --hills HILLS --mintozero --stride 50000 --outfile fes/fes_ttmetad_overall