# OCDesign
![OCDesign method overview](images/OCDesign.png)
OCDesign works as following steps:
* Step 1: Design protein sequences (NeuralMRF, ProDESIGN-LE, or other methods you like)
* Step 2: Assess your proteins for the properties (Solubility, Energy, or all other properties you need)
* Step 3: Select pareto fronts for the next wet-lab validation
* Step 4: Feedback From the wet-lab to correct your Step 1-3


## NeuralMRF: Design protein sequence for a given protein structure
* Input: Protein structure (.pdb file)
* Output：A .fasta file with designed protein sequence

### Environment
Using conda and pip to create the required enviroment for this projection
```
conda env create -f environment.yaml
conda activate neuralmrf
pip install torch==1.10.1+cu111 torchvision==0.11.2+cu111 torchaudio==0.10.1 -f https://download.pytorch.org/whl/cu111/torch_stable.html
```

### Key parameter of the run_neuralmrf.py
```
--seed: A random seed can influence the sequence we design 
--checkpoint_path: model file
--device: which rank gpu to run neuralmrf
--pdb_path: input pdb file
--chain: select chain in the pdb file as the input 
--fix_file: file that specify fixing residues in some positions
--output_file: write the design sequence into which file, Not required
```


### Execute example and result we got
```
# run neuralmrf to design protein sequence that fix residues as we did in round 1
python run_neuralmrf.py --seed 37 --chain A --fix_file example/fix_file/round1.txt --pdb_path example/pdb/5u4y.pdb --device 1
================================================================================
>5u4y Identity:0.4339622641509434
SDNALQNAMKEIQHLPNLDEAEKNSFLLALVLDPSAAEVLRAEARQINIDRQP
================================================================================

# run neuralmrf to design protein sequence that fix residues as we did in round 2
python run_neuralmrf.py --seed 37 --chain A --fix_file example/fix_file/round2.txt --pdb_path example/pdb/5u4y.pdb --device 1
================================================================================
>5u4y Identity:0.660377358490566
FNKAQQNAFYEILHLPNLDEAQKNSFILRLKLDPSAAEVLRAEARQINIDQAP
================================================================================

# run neuralmrf to design protein sequence that fix residues as we did in round 3
python run_neuralmrf.py --seed 37 --chain A --fix_file example/fix_file/round3.txt --pdb_path example/pdb/5u4y.pdb --device 1
================================================================================
>5u4y Identity:0.660377358490566
FAKAQQNAFYEILHLPNLTEEQKNWFILRLKLDPSVAEVLRAEARQINIDQAP
================================================================================
```

## Pareto front selecetion
```
# python pareto_frontier.py CSVFILE 'Property1:[min|max],Property2:[min:max]'  
python  pareto_frontier.py  example/pareto/r2_pareto.csv --sort 'RMSD_VAR:min,MMGBS:min,Solubility:max'  
```
