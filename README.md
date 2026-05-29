# OCDesign: Objective curriculum-guided design of multi-property proteins
![OCDesign method overview](images/OCDesign.png)

## Overview

OCDesign is an objective curriculum-guided computational framework for multi-property protein design. Instead of optimizing all properties simultaneously, OCDesign introduces design objectives sequentially according to a predefined objective curriculum. In each stage, candidate sequences are generated, evaluated across relevant properties, selected using Pareto-front analysis, and then prioritized for experimental validation.

This repository contains code, example data, designed sequences, property evaluation results, and figure-related files used in our study:
### Objective curriculum-guided design of multi-property proteins

| **Stage**     | **Objective(s)**                          | **Purpose**                                          |
|-----------|---------------------------------------|--------------------------------------------------|
| Stage 1   | Solubility and structural consistency | Establish foldable and expressible designs       |
| Stage 2   | Binding affinity                      | Introduce antibody-binding function              |
| Stage 3   | Alkaline resistance                   | Improve robustness under alkaline conditions     |


This objective order was predefined based on domain knowledge to progressively introduce increasingly specific functional requirements.

## OCDesign workflow

**OCDesign consists of four main steps:**

 1. Protein sequence generation
Generate candidate protein sequences using NeuralMRF, ProDESIGN-LE, or other sequence design models.
Property assessment
2. Evaluate designed proteins using computational property predictors or scoring functions, such as solubility, structural consistency, energy, binding-related scores, or task-specific properties.
3. Pareto-front selection
Select candidate designs based on Pareto-optimal trade-offs among multiple properties.
4. Wet-lab validation and curriculum update
Experimentally validate selected candidates and use the results to guide subsequent design decisions and introduce the next objective in the curriculum.

## Installation
We recommend using conda to create the required environment.

```
conda env create -f environment.yaml
conda activate neuralmrf
```
You need install PyTorch (the CUDA 11.1 version) manually
```
pip install torch==1.10.1+cu111 torchvision==0.11.2+cu111 torchaudio==0.10.1 -f https://download.pytorch.org/whl/cu111/torch_stable.html
```

## NeuralMRF: sequence design for a given protein structure

NeuralMRF designs protein sequences conditioned on an input protein structure.

### Input
Protein structure file in .pdb format
Optional fixed-residue file specifying residues that should remain unchanged during design

### Key parameters
```
--seed             Random seed for sequence generation
--checkpoint_path  Path to the trained NeuralMRF model checkpoint
--device           GPU device ID
--pdb_path         Path to the input PDB file
--chain            Chain ID used for sequence design
--fix_file         File specifying fixed residues
--output_file      Output FASTA file; optional
```
### Example: sequence design for each objective stage
**Stage 1: solubility and structural consistency**
```bash
python run_neuralmrf.py \
  --seed 37 \
  --chain A \
  --fix_file example/fix_file/round1.txt \
  --pdb_path example/pdb/5u4y.pdb \
  --device 1
```
Example output:

```bash
>5u4y Identity:0.4339622641509434
SDNALQNAMKEIQHLPNLDEAEKNSFLLALVLDPSAAEVLRAEARQINIDRQP
```


**Stage 2: binding affinity**

```bash
python run_neuralmrf.py \
  --seed 37 \
  --chain A \
  --fix_file example/fix_file/round2.txt \
  --pdb_path example/pdb/5u4y.pdb \
  --device 1
```

Example output:

```bash
>5u4y Identity:0.660377358490566
FNKAQQNAFYEILHLPNLDEAQKNSFILRLKLDPSAAEVLRAEARQINIDQAP
```

**Stage 3: alkaline resistance**

```bash
python run_neuralmrf.py \
  --seed 37 \
  --chain A \
  --fix_file example/fix_file/round3.txt \
  --pdb_path example/pdb/5u4y.pdb \
  --device 1
```

Example output:

```bash
>5u4y Identity:0.660377358490566
FAKAQQNAFYEILHLPNLTEEQKNWFILRLKLDPSVAEVLRAEARQINIDQAP
```
## Pareto-front selection

Candidate designs can be selected based on Pareto-optimal trade-offs among multiple properties.

**Usage**

```bash
python pareto_frontier.py CSVFILE --sort 'Property1:min,Property2:max'
```

Example

```bash
python pareto_frontier.py example/pareto/r2_pareto.csv \
  --sort 'RMSD_VAR:min,MMGBS:min,Solubility:max'
```

In this example:

 - RMSD_VAR:min selects designs with lower structural variation.
 - MMGBS:min selects designs with more favorable binding energy.
 - Solubility:max selects designs with higher predicted solubility.

The script outputs candidate designs on the Pareto front according to the specified property directions.

# Reproducing paper results

The repository provides example files for reproducing the main computational steps:

- Generate sequences using NeuralMRF.
- Assess designed sequences using property scores provided in data/property_scores/.
- Perform Pareto-front selection using pareto_frontier.py.
- Compare selected designs with experimentally validated candidates.
- Reproduce figure panels using files in results/figures/.

## Citation

If you use OCDesign or the data in this repository, please cite:

Objective curriculum-guided design of multi-property proteins.
Manuscript under review.

The citation will be updated after publication.


## Contributors
- Longying Liu., Jianquan Zhao, and Xiaomin Xie  contributed equally.
- Xinmiao Liang, Xianlong Ye, Dongbo Bu., and Han Zhou. conceptualized the study.
- Longying Liu. performed wet lab experiments, data curation, and wrote the original draft.
- Jianquan Zhao. performed multi-round molecular design.
- Xiaomin Xie  performed molecular simulation design and MD simulation-based screening.
- Simeng Xu. assisted with wet lab experiments.
- Xinru Zhang. performed primary protein sequence design using ProDESIGN-LE.
- Zaikai He. revised the manuscript.
- Chungong Yu. provided top-level design and guidance.
- Fan Liu. assisted with molecular simulation.
- Kun Wang. contributed to protein sequence design using ProDESIGN-LE.
- Xinglong Wang. contributed to  NeuralMRF development.
- Milong Ren. revised the manuscript.
- Xinmiao Liang. provided top-level design and guidance.
- Xianlong Ye. supervised the study, provided overall conceptualization and guidance.
- Dongbo Bu. supervised the study, provided overall conceptualization and guidance, and revised the manuscript.
- Han Zhou. supervised the study, provided overall conceptualization and guidance, and revised the manuscript.

## License
- Code: MIT License or Apache License 2.0
- Data and figures: CC BY 4.0,
## Contact
For questions, please contact:

Dongbo Bu (Corresponding author): dbu@ict.ac.cn
Jianquan Zhao (Co-first author): zhaojianquan@ict.ac.cn

