# HeteroGRN: Gene Regulatory Network Inference with Heterogeneous Graph Neural Networks

## Requirements
Dependencies are provided in the following files:
- Python 3.13.2: requirements.txt
- R 4.3.3: requirements_R.txt

## Usage

### 1. Main Data Source
- BEELINE benchmark dataset (scRNA-seq data and ground-truth networks) is available at:
  https://zenodo.org/records/3701939
- Protein structure data (predicted 3D structures) is downloaded from AlphaFold DB:
  https://alphafold.ebi.ac.uk/
- For missing protein structures, supplementary structures can be obtained from UniProtKB:
  https://www.uniprot.org/uniprotkb

### 2. Code Structure
- ./code/datapreprocess/: Data preprocessing scripts, used to generate training/validation/test sets and multi-source similarity matrices (expression, GO, protein structure).
- ./code/main.py: Main script for model training and testing.
- ./Demo/demo.py: Script for real-world application and GRN prediction.

### 3. Important Notes
- All file paths in the code should be modified according to your actual data storage location.
- Default parameters are recommended for low-density networks.
- For high-density networks, use the following parameters:
--dropout_rate 0.45 \
--edge_dropout_rate 0.55 \
--gnn_hidden 512 \
--gnn_out 32 \
--patience 40
### 4. Evaluation (Training & Testing)
Run the main script to train and evaluate the model on the BEELINE dataset.
Example Command (Low-density network)
python ./code/main.py \
  --cell_type hESC \
  --net_type NonSpecificDataset \
  --gene_num 500 \
  --fold 1
Example Command (High-density network)
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
Outputs will be saved to:
·./example/model/ : trained model files
·./example/val/ : validation results
·./example/test/ : test results 
### 5. Application (Demo Prediction)
We present an application example of the model for predicting gene regulatory relationships.

Prepare Data
All required files for the application example are provided in the ./Demo/ directory, including training set, validation set (for early stopping), gene pairs to predict, and multi-source similarity matrices.

Run the Demo
cd Demo
python demo.py

Output
Predicted GRN results will be saved to:
./Demo/application_results/
