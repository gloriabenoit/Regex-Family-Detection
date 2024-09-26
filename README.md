# Protein family pattern detection using regular expressions
December 2022 (L3 BI)

## Introduction

This project aims to detect a protein's family from a signature found in its sequence, using regular expressions.
All protein domains and families are found on the [PROSITE](https://prosite.expasy.org/) database. 
PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.

Starting from a *PROSITE* motif, we have automated its translation into a regular expression so that it can be searched in sequences. 

## Setup

To install the algorithm and its dependencies, you need to perform the following steps:

### Clone the repository

```bash
git clone https://github.com/gloriabenoit/Regex-Family-Domain.git

cd Regex-Family-Domain
```

### Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Create a Conda environment

```bash
conda env create -f environment.yml
```

### Activate the Conda environment

```bash
conda activate regex-family-detection
```

## Usage

`Detection_MotifProtéine.ipynb` takes you through everything you need to know about the project, while `demo.py` is pure code.

