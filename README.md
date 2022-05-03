# Transgenerational DNA Damage Code
The repository contains the custom code used in the paper Transgenerational inheritance of paternal DNA damage by linker histone H1-mediated DNA repair restriction. To detect de novo genomic aberrations we performed WGS of single P0 and their direct F1 progeny upon IR irradiation. Structural variants were called with Manta v1.6.0 and the resulting VCF files are used as input for this analysis. 
preprocess.py does some quality filtering and marks F1 SVs that already appeared in the P0 generation. 
analysis.py contains the code for the relevant figures in the paper.
