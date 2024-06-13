# CB2M_scTransformers 
**Project to analyze the behavior and prediction quality of single-cell annotation tools based on Transformers** 

## Table of Contents 
1. [Data](#data)
2. [scBERT](#scbert)
   - [Installation](#installation)
   - [Usage](#usage)
3. [scGPT](#scgpt)
   - [Usage](#usage-1)
4. [To-Do List](#to-do-list)

## Data 
Detailed data description and preprocessing steps will be provided soon.

## scBERT 
![Python 3.8+](https://img.shields.io/badge/python-3.6.8-brightgreen) 
### Installation 
![scipy 1.5.4](https://img.shields.io/badge/scipy-1.5.4-yellowgreen) ![torch 1.8.1](https://img.shields.io/badge/torch-1.8.1-orange) ![numpy 1.22](https://img.shields.io/badge/numpy-1.19.2-red) ![pandas 1.1.5](https://img.shields.io/badge/pandas-1.1.5-lightgrey) ![scanpy 1.7.2](https://img.shields.io/badge/scanpy-1.7.2-blue) ![scikit-learn 0.24.2](https://img.shields.io/badge/scikit__learn-0.24.2-green) ![transformers 4.6.1](https://img.shields.io/badge/transformers-4.6.1-yellow) ![matplotlib 3.6](https://img.shields.io/badge/matplotlib-3.6-blue) 
### Usage 

#### 01_Datapreprocessing 
The first preprocessing phase applies to both scBERT and scGPT. Its goal is to sort, rename, and select values to be retained for use by these tools. 

**Output:** 
- 1 file H5AD (Example : `X_Preprocess.h5ad`)
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`
- `metadata.csv`

#### 02a_GlobalHeterogenity (Variable Gene)

To obtain the varabiable genes from R, here's the procedure: apply the entire code in the following order:
- `analysisParams`
- `00_generalDeps`
- `01_prepareData`
- `02_variableGenes`

**Output:** 
- `Variable_Gene.csv`

#### 02b_FilterData 
Unlike scGPT, scBERT requires Bash for execution as it calls Python functions from Bash. Preprocessed files for thet training and the best of cell types must be output beforehand. 

**Output:** 
- Double the number of selected folds, providing test and training datasets (e.g., if 5 folds are chosen, there will be 5 test datasets and 5 training datasets).
- Example : `X_test_FOLD_0_Preprocess.h5ad` `X_Training_FOLD_4_Preprocess.h5ad`

Note that the file named "training" is to be used for the fine-tuning phase, and the "test" file is to be used for predictions.

#### 06_scBERT
A Jupyter Notebook for running scBERT predictions using the preprocessed data. 

**Output:** 
- 3 directories (epoch 1, 5, and 10) containing datasets with:
- Example : `Cell_Type_Predict_FOLD_0_epoch_1` = File containing the cell type predict
- `label` and `label_dict` = File was used for the fine-tuning and create per scBERT himself.
- `finetune_best.pth` = File model create, and use for the prediction of the cell types.

#### 07_scBERT_Analysis 
A Jupyter Notebook for analyzing the prediction results, including indicator calculations. 

## scGPT 
![Python 3.9+](https://img.shields.io/badge/python-3.6.8-brightgreen) 

### Installation 
![torch 1.8.1](https://img.shields.io/badge/torch-1.8.1-orange) ![markupsafe 2.0.1](https://img.shields.io/badge/numpy-1.19.2-red) ![scgpt](https://img.shields.io/badge/pandas-1.1.5-lightgrey) ![wandb](https://img.shields.io/badge/scanpy-1.7.2-blue) ![packing](https://img.shields.io/badge/scikit__learn-0.24.2-green) ![transformers 4.6.1](https://img.shields.io/badge/transformers-4.6.1-yellow) ![matplotlib 3.6](https://img.shields.io/badge/matplotlib-3.6-blue) [scikit-learn 0.24.2](https://img.shields.io/badge/scikit__learn-0.24.2-green)

### Usage 

#### 01_Datapreprocessing 

If the phase has not been previously performed, then follow the same procedure described above. 

Note: If you've already done this step with scBERT or scGPT on the same dataset, you don't need to do it again. You can re-use the previously obtained files when they are needed in the code later.

#### 02a_GlobalHeterogenity (Variable Gene)

If the phase has not been previously performed, then follow the same procedure described above. 

Note: If you've already done this step with scBERT or scGPT on the same dataset, you don't need to do it again. You can re-use the previously obtained files when they are needed in the code later.

#### 03_scGPTAnalysis

Unlike scBERT, folds and dataset selection and preprocessing are done directly in the Jupyter Notebook. 

**Output:** 

- `indicators_results_annotation`: Indicator values calculated for all folds (blocks).
- `prediction_result_annotation`: File containing the IDs of the cell types predicted for each cell.
- `id_map_results_annotation`: File containing cell type IDs.
- `dataframe_liste_number_per_fold`: File containing the number of cells present in folds by cell type.
- `labels_result_annotation`: File containing the actual cell types for each cell in the dataset.
- `confusion_matrix` : png and csv.


For analysis, I decided to split it into two parts, which can be combined (as in the case of scBERT).

#### 04_scGPTAnalysis_Result

Analysis of results from files containing a single epoch (one run). 

**Output:** 
- HTML file by Quarto 

#### 05_scGPT_Different_Epoch 

Analysis of results across different epochs (runs). 

**Output:** 
- HTML file by Quarto

## To-Do List 
- [ ] Automate the output directories for scBERT
- [ ] Obtain HTML outputs via Quarto for scBERT 
- [x] Test 
