#!/bin/bash

prefix=$1

matlab -nodisplay -nodesktop -r "cd('/cbica/projects/pinesParcels/multiscale/scripts/derive_parcels/Toolbox/PBP/'); PBP_vertWiseEffect_Bart('/cbica/projects/alpraz_EI/output/surface_files/${prefix}_LH.csv','/cbica/projects/alpraz_EI/output/surface_files/${prefix}_RH.csv','/cbica/projects/alpraz_EI/output/surface_files/${prefix}.svg');exit"