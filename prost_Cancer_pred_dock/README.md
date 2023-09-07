# prost_Cancer_pred_dock
This is a simple platform for computing Lipinsky's Rule of five using the rdkit package, predicting pIC50 of canonical SMILES that are potential targets against Oestrogen receptor alpha protein as ant-prostate cancer agaents using apreformatted RandomForest model, and docking of the canonical SMILE with the Oestrogen receptor alpha protein using Audodock Vina package. 

### Purpose of the Package
+ The purpose of the package is to provide a unified platform for computing prostate cancer drug likeness indicess and performing docking on the same compounds. 

### Features
+ Important chemoinformatics features of Oestrogen receptor alpha antagonists such as:
    - Lipinsky descriptors
    - Prediction of pIC50
    - Docking and visiualization 

### Getting Started
The package is found on pypi hence can be installed with pip

#### Installation
```bash
pip install prost_cancer_pred_dock
```
#### Dependencies
+ pymol, rdkit, pandas, padelpy, joblib, csv, meeko, Autodock Vina, java, scipy, numpy, and scikit-learn.
+ However, most of them can be easily installed with conda and pip.     
```bash
conda install -c conda-forge pymol-open-source
conda install -c bioconda autodock-vina
conda install -c cyclus java-jre
pip install rdkit-pypi pandas padelpy joblib python-csv meeko scikit-learn scipy numpy
```

### Usage
#### Computation of Lipinsky descriptors
```python
from prot_cancer_pred_dock import calculate_lipinski_descriptors
calculate_lipinski_descriptors("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
```
#### Prediction pIC50
```python
from prost_cancer_pred_dock import predict_pIC50
predict_pIC50("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
```
#### Docking and visualization
```python
from prost_cancer_pred_dock import prot_lig_docking
prot_lig_docking("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
```
#### Visualization of docking results
```python
from prost_cancer_pred_dock import vizualize_dock_results
vizualize_dock_results()
```


### Contribution
We welcome any contributions. Should you notice a bug, please let us know through issues in the github repository. Thank you very much.

### Authors
+ Edwin mwakio
+ Clabe Wekesa
+ Patrick Okoth
