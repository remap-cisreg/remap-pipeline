# remap-pipeline

## Description
ReMap is a project which goal is to provide the largest catalogue of high-quality regulatory regions resulting from a large-scale integrative analysis of hundreds of transcription factors and general components of the transcriptional machinery from DNA-binding experiments.
This git contain all files necessary to run analysis ReMap style from annotaion files to the final product !


## Installation
### Requirements
 - Python 3
 - Snakemake >= 5.5.1
 - Recommended : Conda/Docker/Singularity
### Step by step
 1. Pull the git
 2. That's all !
 
## Usage
More in depth info are in the [wiki](https://github.com/remap-cisreg/remap-pipeline/wiki).

### General

 1. repare the metadata from your annotation file by extracting downloading info (See README.md in 1.metadata/)
 2. Create cluster config and snakemake config to match your set up (See example in 2.scripts/cluster_configuration)
 3. Get necessary files such as reference genome (See README.md in 3.genome/)
 4. Create a launch bash script (See example in root directory)
 5. Run launch script
 
### Compute clusters specificities
remap-pipeline can be run with Conda or Docker/Singularity and Torque and Slurm.
Examples are in launch scipts.

## Contributing

- Jeanne Cheneby
- Lionel Spinelli
- Benoit Ballester

## License
Under GNU GPLv3 licence.

## References
Manuscrit in preparation
