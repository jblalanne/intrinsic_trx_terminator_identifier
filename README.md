# intrinsic_trx_terminator_identifier

These scripts identify putative intrinsic terminators by looking for long stretches of Ts downstream (corresponding to U-tracts) of genes preceeded by strong and clean RNA hairpins. 

The pipeline takes in .fna (genome sequence) and .gff (gene annotation) files (examples for Bacillus pseudofirmus included in example_GFF_fasta_files), and output positions and stop-to-stem distances for putative terminators. 

RNAfold must be installed. The Matlab script calls RNAfold and reads output directly. Other dependencies are in the subroutines directory. The directory path to batch_RNAfold.txt in line 18 of master_RNAfold_file_v2.m will also need to be updated.

To run the script, first run terminator_identification_header_script_v4_git_20200504.m and terminator_identification_header_random_sequence_git_20200504.m. These respectively identify putative RNA hairpins upstream of Ts close to ORFs and fold random positions in the genome. To threshold on hairpins parameters using the output of the two above scripts, parsing_hairpin_data_structures_v2_git_20200504.m must be run. To parse the final data structure, generate_final_term_data_structure_git_20200504.m can be used. Example of data variable output for Bacillus pseudofirmus can be found in example_output_variables). 
