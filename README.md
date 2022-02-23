# RDF
Calculate radial distribution function with C++ backend and python frontend
https://github.com/ankur9511/fortan_dcd_analysis

Use above to convert xtc to dcd
module load python/2.7.15

module load gcc

module load intel

git clone https://github.com/ankur9511/fortan_dcd_analysis.git

cd fortan_dcd_analysis/

chmod +x ./script

./script "createmodule"

python calculations.py
