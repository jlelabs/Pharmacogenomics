# Pharmacogenomics

This repository outlines an end-to-end pipeline using Stargazer, bcftools and an in-house Python script to determine metabolizer phenotypes from PharmGKB and CPIC based on genotyping. 

## Built With
- Stargazer
- bcftools
- Python

## Getting Started
### Installation
```
git clone https://github.com/jlelabs/Pharmacogenomics
```
### Configuration Files
Download the following from Stargazer https://stargazer.gs.washington.edu/stargazerweb/res/form.html and transfer them into your Configuration folder within the Pharmacogenomics repository. 
- Stargazer_v1.0.8

### Usage

Download and configure NA12878 sample
```
./NA12878.sh
```

Run main Pharmacogenomics script
```
./Pharmacogenomics.sh -i ${SAMPLEID} -v ${VCF} -b ${BAM} -o ${OUTPUT_DIRECTORY}
```
