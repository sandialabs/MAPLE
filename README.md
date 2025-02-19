# MAPLE

##Note: this repository was formerly named EAGLE

## Usage
Implements a tensorflow mixed CNN and multi-head attention model for predicting gene expression from epigenetic data within or across organism species. The codebase includes data preparation helper functions (set_TSS_window_index.py and crosscompcov.py) that take pre-processed .bam files as inputs and provides data transformed for use by MAPLE for use with new data. Note that the code does set a seed for reproducibility. 

## Environment Setup Dependencies
The recommended setup is to clone this repository and to use conda to create the environment that MAPLE will run in. Use:
`conda env create -f MAPLE-env.yml`
and 
`conda activate MAPLE-env`

## For Use With Data from Publication
Those interested in using the species and data already included in MAPLE may start directly with `MAPLE_main.py`, which includes its data loader.

Positional Arguments:
- `species1`: "Ncrassa", "Fgram", "Anid", or "LmacL". Default is "Ncrassa"
- `species2`: "Ncrassa", "Fgram", "Anid", or "LmacL". Default is "Fgram". It is possible to run without a second species by including an empty string ("").
- `n_outputs`: Number of bins for predicted gene expression classification. Default is 2.
- `bin_num`: Size of kernel that input data will be binned into. Smaller bin_num will have more resolution 
of input data signal profile. Default is 50.

Keyword Arguments

- `--epochs`: Number of epochs to train with. Default is 100.
- `--batch_size`: Batch size for training. Default is 100.
- `--trained_model_result`: If provided, the trained model will be saved to the provided "path/file" (type=str).
- `--runHPO`: If True, a grid search hyperparameter optimization will be run to train and optimize the model. Default is False.

Using 
`model = MAPLE.py(args)`
will return the trained model object for provided species, which may be used for evaluation tasks.

## For Use With Custom Data
It is possible to prepare and use new data with MAPLE. To do so, the developer will need to run pre-processed (e.g., QC checked, trimmed, aligned, indexed, and sorted) bam files through the helper functions before using MAPLE.
Order of usage:
1. RNA_norm (python function)
2. set_TSS_window_index.py (command line)
3. crosscompcov.py (command line) 
4. MAPLE_main.py, which uses loaddeepdata2.py (command line)

### Usage of RNA_norm

    from Tmm_normalization import RNA_norm
    RNA_norm(rawRNAdata, normalizedRNAdata)
With input arguments:
- `rawRNAdata`: path to file preprocessed RNA-seq data
- `normalizedRNAdata`: filename to write normalized and averaged RNA-seq data to

This will write the normalized RNA-seq data to the filename and path provided in the normalizedRNAdata argument. It is recommended to condition match the RNA-seq samples used with the epigenetic data being used, and to pass only those columns into this function.

### Usage of set_TSS_window_index.py
Processes gene annotations to define windows around transcription start sites. If
input(s) are GFF3, convert to BED.

Using the gff2bed library, convert data from GFF3 to BED. List input file(s) as
arguments, and the outputs will be written to files with the BED extension:

    python set_TSS_window_index.py alpha.gff subdir/bravo.ext charlie.bed
will write: 

    alpha.bed
    alpha.promoter.bed
    subdir/bravo.promoter.bed 
    charlie.promoter.bed 

### Usage of crosscompcov.py
Gets epigenetic coverage upstream and downstream of a transcription start site, using bed file from set_TSS_window_index.py.

To use, pass in a BED file with features to check, and a collection of read depth files
to check for read coverage statistics.

    python crosscompcov.py --bed promoters.bed organism/*.depths.txt > promoter.coverage.p

### Usage of MAPLE_main.py
In addition to the arguments defined above, providing the following keyword arguments will allow the user to pass in prepared data for any organism.

Keyword Arguments:
- `--species1cov`: relative path to file with epigenetic coverage in peri-transcription site window for species 1
- `--species1rna`: relative path to file with RNA expression data for species 1
- `--species1mods`: list of strings defining epigenetic modifications for species 1 
- `--species2cov`: relative path to file with epigenetic coverage in peri-transcription site window for species 2
- `--species2rna`: relative path to file with RNA expression data for species 1
- `--species2mods`: list of strings defining epigenetic modifications for species 2

These will all default to None if not provided. If providing one custom species dataset, both species must be provided. If only  species names are used, the first positional arguments, it will use the data from the original publication.

The user must provide pickle files output from crosscompcov.py for each species coverage file and the RNA reads files as .txt files, ideally normalized. The coverage and RNA read files must be annotated from the same reference genome so the genes are consistently named. Further, the user must specify the epigenetic modifications to be used for each species, which must match the column headers in the coverage .p file and they must be shared between the species being called, in the same order. Paths must be relative to working directory. The coverage files must have gene names in a column named "Gene", and RNA data need to be in an averaged column called "Average". Use of the RNA_norm function provided here will produce a file suitable for use with MAPLE.

## Planned Improvements
- Increased support for using custom data
- Extraction of gene names from data loader

## Citing MAPLE
Once a preprint is available, the citation will be here.
