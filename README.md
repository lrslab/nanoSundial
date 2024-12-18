# nanoSundial
A de novo modification detection tool targeting current features based on the 004kit for prokaryotes.

## Environment preparation

In order to make it easy to run nanoSundial, we provided two different methods for users to prepare the software environment.

(1). Installing packages from **Conda**

    conda install samtools=1.16 minimap2 f5c=1.4 slow5tools -c conda-forge -c bioconda 
    git clone https://github.com/lrslab/nanoSundial
    cd nanoSundial/
    pip install -r requirement.txt

(2). Installing from **Docker**

    docker pull zhihaguo/nanocem_env
    # git clone is required when it enters to the Docker container

## Pre-processing

In theory, nanoSundial is a comparative method that requires **two** samples. 

Using a **wild type** and a **negative control**is recommended. 
Raw sequencing data (**blow5**) and basecalled file (**fastq**) files are needed for each sample. 
For the content of how raw data convert to blow5 and how to basecall, please refer to [this](https://nanocem.readthedocs.io/en/latest/preparation/)

Additionally, a reference sequence (**fasta** file) is required. The genome is recommended for prokaryotic analysis; however, the transcriptome might be more suitable for eukaryotic analysis.
The file list should be structured as follows, and q quick start with example data is [here](#quick-start-with-example-data):

    test_data/
    ├── wt/
    │   └── sample.blow5  # Raw signal, BLOW5 file  
    │   └── sample.fastq # Basecalled reads, FASTQ file
    ├── ivt/
    │   └── sample.blow5  # Raw signal, BLOW5 file   
    │   └── sample.fastq # Basecalled reads, FASTQ file
    └── ref.fasta  # Reference sequence, FASTA file

## Manual of nanoSundial
nanoSundial offers three functional scripts: `extract_feature_region.py`, `sundial_comp.py`, and `merge_positive_region.py`. 
The extract_feature_region.py script needs to be run separately for each of the two samples, ensuring the same parameters are used for each run.
### 1. Feature extraction
We have provided a script for feature extraction called `extract_feature_region.py`. 

Below is an introduction to the parameters

    python extract_feature_region.py -h
    usage: extract_feature_region.py [-h] [-f FASTQ] [-b BLOW5] [-o OUTPUT]
                                     [-r REF] [-t CPU] [--pore {r9,r10,rna004}]
                                     [--subsample SUBSAMPLE] [--win_size WIN_SIZE]
                                     [--flip FLIP] [--rna]
                                     [--base_shift {auto,0,-1,-2,-3,-4,-5,-6,-7,-8}]
    
    optional arguments:
      -h, --help            show this help message and exit
      -f FASTQ, --fastq FASTQ
                            input fastq file
      -b BLOW5, --blow5 BLOW5
                            input blow5 file
      -o OUTPUT, --output OUTPUT
                            output file path (will be created if this folder does not exist, default:nanosundial_feature)
      -r REF, --ref REF     reference path (fasta file)
      -t CPU, --cpu CPU     Process numbers
      --pore {r9,r10,rna004}
                            pore (default: rna004)
      --subsample SUBSAMPLE
                            subsample ratio (0-1, used to subsample data, default:1)
      --win_size WIN_SIZE   windows size (default:100000)
      --flip FLIP           length of the base of exclusion of start or end
      --rna                 Turn on the RNA mode
      --base_shift {auto,0,-1,-2,-3,-4,-5,-6,-7,-8}
                            base shift option (default: auto)

### 2. Statistical analysis
After feature extraction from two samples, `sundial_comp.py` is provided to compare the two samples. 
We offer three algorithms: Multivariate Analysis of Variance (**manova**), Logistic Regression (**lr**), and the Kolmogorov-Smirnov test (**ks**), with `MANOVA` as the default. 
The script returns each position that meets the coverage criteria (Default:**50x**).

Below is an introduction to the parameters

    python sundial_comp.py -h
    usage: sundial_comp.py [-h] [-i INPUT] [-c CONTROL] [-o OUTPUT] [-t CPU] [--balance] [--method METHOD]
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            output folder of WT from extract_feature_region.py
      -c CONTROL, --control CONTROL
                            output folder of negative control from extract_feature_region.py
      -o OUTPUT, --output OUTPUT
                            output file (default:nanosundial_result.csv)
      -t CPU, --cpu CPU     Process numbers
      --balance             Turn on the balance mode
      --method METHOD {'manova', 'lr', 'ks'}
                            statistical algorithms provided(default:manova)
      --minimum_coverage MINIMUM_COVERAGE
                            coverage cutoff of this comparison (default:50)


### 3. Merging strategy
After applying the recommended cutoff values (|mean difference|: 0.18, |dwell time difference|: 1, -log10(adj p-value): 3), all points within a distance less than shift_size are connected into a region via `merge_positive_region.py`. 

    python merge_positive_region.py -h
    usage: merge_positive_region.py [-h] --input INPUT --output OUTPUT [--coverage_cutoff COVERAGE_CUTOFF] [--mean_cutoff MEAN_CUTOFF] [--dwell_cutoff DWELL_CUTOFF]
                         [--pvalue_cutoff PVALUE_CUTOFF] [--shift_size SHIFT_SIZE]
    
    Merging the raw result of nanoSundial
    
    optional arguments:
      -h, --help            show this help message and exit
      --input INPUT         Input CSV file
      --output OUTPUT       Output BED file (default:nanoSundial_merged_region.bed)
      --coverage_cutoff COVERAGE_CUTOFF
                            Coverage cutoff value (default:50)
      --mean_cutoff MEAN_CUTOFF
                            cutoff of absolute mean difference (default:0.18)
      --dwell_cutoff DWELL_CUTOFF
                            cutoff of absolute dwell time difference (default:1)
      --pvalue_cutoff PVALUE_CUTOFF
                            cutoff of -log10(adj p-value)  (default:3)
      --shift_size SHIFT_SIZE
                            shift size of merging adjacent position (default:4)

## Quick start with example data

    # Download data
    git clone https://github.com/lrslab/nanoSundial.git
    cd nanoSundial/
    # Feature extraction
    cd test_data/wt/
    python ../../extract_feature_region.py --fastq sample.fastq --blow5 sample.blow5 --ref ../ref.fasta --rna --cpu 4
    cd ../ivt/
    python ../../extract_feature_region.py --fastq sample.fastq --blow5 sample.blow5 --ref ../ref.fasta --rna --cpu 4

    # Statistical analysis
    cd ../
    python ../sundial_comp.py -i wt/nanoSundial_feature/ -c ivt/nanoSundial_feature/ -t 4

    # Merging strategy
    python ../merge_positive_region.py --input nanosundial_result.csv
