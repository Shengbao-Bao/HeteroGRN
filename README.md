# HeteroGRN: Gene Regulatory Network Inference with Heterogeneous Graph Neural Networks

## Requirements

Dependencies are listed in the following files:

- **Python 3.13.2**: `requirements.txt`
- **R 4.3.3**: `requirements_R.txt`

---

## Usage

### 1. Main Data Sources

- **BEELINE benchmark dataset** (scRNA-seq data and ground-truth networks):  
  https://zenodo.org/records/3701939

- **Protein structure data** (predicted 3D structures):  
  AlphaFold DB  
  https://alphafold.ebi.ac.uk/

- **Supplementary protein structures for missing entries**:  
  UniProtKB  
  https://www.uniprot.org/uniprotkb

---

### 2. Code Structure

- `./code/datapreprocess/`  
  Data preprocessing scripts used to generate training, validation, and test sets, as well as multi-source similarity matrices (expression, GO, and protein structure).

- `./code/main.py`  
  Main script for model training and evaluation.

- `./Demo/demo.py`  
  Script for real-world application and GRN prediction.

---

### 3. Important Notes

- Please modify all file paths in the code according to your local data storage locations.
- Default parameters are recommended for **low-density networks**.
- For **high-density networks**, use the following parameters:

```bash
--dropout_rate 0.45 \
--edge_dropout_rate 0.55 \
--gnn_hidden 512 \
--gnn_out 32 \
--patience 40
```

---

### 4. Evaluation (Training and Testing)

Run the main script to train and evaluate the model on the BEELINE dataset.

#### Example command for a low-density network

```bash
python ./code/main.py \
  --cell_type hESC \
  --net_type NonSpecificDataset \
  --gene_num 500 \
  --fold 1
```

#### Example command for a high-density network

```bash
python ./code/main.py \
  --cell_type mDC \
  --net_type SpecificDataset \
  --gene_num 500 \
  --fold 1 \
  --dropout_rate 0.45 \
  --edge_dropout_rate 0.55 \
  --gnn_hidden 512 \
  --gnn_out 32 \
  --patience 40
```

#### Output directories

Outputs will be saved to:

- `./example/model/` — trained model files
- `./example/val/` — validation results
- `./example/test/` — test results

---

### 5. Application (Demo Prediction)

We provide an application example for predicting gene regulatory relationships using our model.

#### Prepare Data

All required files for the application example are provided in the `./Demo/` directory, including:

- training set
- validation set (for early stopping)
- gene pairs to predict
- multi-source similarity matrices

#### Run the Demo

```bash
cd Demo
python demo.py
```

#### Output

Predicted GRN results will be saved to:

- `./Demo/application_results/`
