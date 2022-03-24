# intrinsic_trx_terminator_identifier

These scripts identify putative intrinsic terminators by looking for long stretches of Ts downstream (corresponding to U-tracts) of genes preceeded by strong and clean RNA hairpins. 

These scripts were used to generate data for [Johnson and Lalanne et al. 2020, Nature](https://www.nature.com/articles/s41586-020-2638-5).

For more information or if you have any questions, post an issue to this repository or get in touch with the [Li Lab](http://gwli.scripts.mit.edu/group/). 

## Input and Output

Input: genome sequence(s) and annotation(s)

Output: positions and stop-to-stem distances for putative terminators.

## Prerequisities 

* A genome sequence (.fna)
* A genome annotation (.gff)

(See `example_GFF_fasta_files`)

* RNAfold (from the ViennaRNA package).
* 
By default, the `subroutines/batch_RNAfold.sh` script calls `usr/local/bin/RNAfold` (line 15) to generate RNAfold. If your path to RNAfold differs, you will need to modify that line.

* Various helper functions (included in `/subroutines/` in this repository, just make sure to add them to your MATLAB path)

* MATLAB + Signal Processing Toolbox + Statistics and Machine Learning Toolbox

The `findpeaks` command requires the signal processing toolbox.

The `hist3` command requires the statistics and machine learning toolbox.

## Running

First, put your .fna genome sequence(s) in the `genomes_fasta` folder. 

Then, put your .gff genome annotation(s) in the `GFFs` folder.

Then, run `run_pipeline.m` (or type `run_pipeline` in the MATLAB command line).

If all goes well, your final results files will appear in a folder called `final_stop_to_stem_results`.

This takes on the order of a few minutes for each bacterial genome. 

## Notes

* On an M1 Macbook with MATLAB 2022a, this pipeline was freezing for me on a plotting step. Downgrading to MATLAB 2019a fixed the problem for me. 

### Pipeline explanation

`terminator_identification_header_random_sequence_git` runs first, which takes random sequences from the genome and runs those through RNAfold to generate a baseline for the number of hairpins you would see by chance.

Then, `terminator_identification_header_script_v4_git` takes sequences downstream of T-rich regions (U-tracts) and runs RNAfold to generate hairpin characteristics.

The amount of hairpins you would see by chance is dependent on GC content, so `get_genome_annotations_GC_content_git` computs genome-wide characteristics to control by.

Finally, `parsing_hairpin_data_structures` and `generate_final_term_data_structure` parses the previously-generated MATLAB variables into an easy-to-export format. 
