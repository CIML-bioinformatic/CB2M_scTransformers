# CB2M_scTransformers 

The annotation of cell types on single cell RNA-seq data is a complex, uncertain and time-consuming task, requiring several methods and tools to be able to annotate cells appropriately and efficiently. To overcome these problems and uncertainties, numerous tools and scientific articles have emerged over the years.
The rise of artificial intelligence in our lives (notably through chatGPT), has also imposed itself on the scientific world, bringing novelty and innovation to existing techniques and tools. These tools need to be tested and studied to verify their effectiveness. In this project, two cell annotation tools in single cell RNA-seq named scBERT and scGPT are of interest to CB2M because of their ability to resolve and avoid the uncertainties and problems mentioned above.
We study here, through various analyses, including cross-validation and the use of multiple qualitative and numerical indicators, that cell annotation by those tools are effective for annotating cells from scRNA-seq.

## Table of Contents 
1. [Goal of the github](#goal-of-the-github)
2. [Description of the datasets](#description-of-the-datasets)
3. [Prepare the environments](#prepare-the-environments)
   - [Clone the github repository](#clone-the-github-repository)
   - [Set the WORKING_DIR variable](#set-the-working_dir-variable)
   - [Add you working dir in the code](#add-you-working-dir-in-the-code)
   - [Download the data](#download-the-data)
   - [Install Singularity and Docker](#install-singularity-and-docker)
   - [Download the Singularity images](#download-the-singularity-images)
   - [Download the Docker images and load them on your system](#download-the-docker-images-and-load-them-on-your-system)
4. [Run the analysis](#run-the-analysis)
   - [Run the analysis workflow](#run-the-analysis-workflow)
   - [Run the analysis individually using Docker](#run-the-analysis-individually-using-docker)
5. [Tools Transformers](#tools-transformers)
   - [scBERT](#scbert)
   - [scGPT](#scgpt)
---
---
## Goal of the github
This github project contains the instructions and material to reproduce the analyses reported in this project.
Source code (scripts and dockerfiles) are available in the github repository. Required data and built Docker/Singularity images are available on download. Instructions to reproduce the analyses are provided below.

To reproduce the analysis, you have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

---

## Description of the datasets

There are 2 datasets in this study. 

* Cell atlas of human thymic development : 15 embryonic and fetal thymuses covering stages of thymic development from 7 post-conceptional weeks (PCW) to 17 PCW, and 9 postnatal thymuses from human pediatric and adult samples. It contains 255,901 cells, 32,922 genes and 33 different cell types.
* Cell atlas across tissue in human immune system : Immune compartment of 15 tissues from six deceased adult donors. It contains 329,762 cells, 36,398 genes and 35 different cell types.

When downloading the code and data (see below), you will obtains 2 sub-folders with names as below:

```
    scTransformers
    │
    ├── cross_tissue_immune_cell
    │
    └── Human_Thymus_Development_Atlas
```
---
---
## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the github repository and set the WORKING_DIR environment variable
- Download the pre-processed data
- Install Singularity and Docker
- Download the Singularity images
- Download the Docker images and load them on your system

Below you will find detailed instruction for each of these steps.

---

### Clone the github repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder **scTransformers** with all the source code. 

---

### Set the WORKING_DIR variable

Then, you must set an environment variable called **WORKING_DIR** with a value set to the path to this folder.

For instance, if you have chosen to clone the Git repository in __"/home/thyarion/workspace"__, then the **WORKING_DIR** variable will be set to __"/home/thyarion/workspace/scTransformers"__.

**On linux:**

```
    export WORKING_DIR=/home/thyarion/workspace/scTransformers
```

---

### Add you working dir in the code

The code uses variables that are stored in different "parameters" file. One important variable is the PATH_PROJECT which indicate to the code where your project is stored.
You have to modify this variable in the code to reflect your project setup. For the step of variable gene, u will need to modify the working dir in **AnalysisParams.R** in the subfolder **03_Script**.

For python, edit in each files the line defining the **PATH_PROJECT** variable and change its value to the same value as the **WORKING_DIR** variable you defined before. Then save the files.

```
PATH_PROJECT = "/home/thyarion/workspace/scTransformers"
```

---

### Download the data

Each sample needs its own sub-folder containing the initial data used by the analysis and the output files of the analysis. Those data can be downloaded from Zenodo and uncompressed. The Zenodo dataset DOI are 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11658106.svg)](https://doi.org/10.5281/zenodo.11658106) (Human Thymus Development Atlas) and [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11658091.svg)](https://doi.org/10.5281/zenodo.11658091) (Cross tissue immune cell)

List of directories in reference files:
* 00_Dataset : contains the H5AD files dowloaded, and the metadata file.
* 01_CellCycle : Files use by R for the selection of variable gene
* 02_HeatShock : Files use by R for the selection of variable gene
* 04_Model : contains the pre trained model of scBERT and scGPT. But also, the files .py use by scBERT to do the prediction.

List of directories in output files :
* 01_Datapreprocessing : contains the result of pre processing data use for the anlysis of variable gene.
* 02a_GlobalHeterogenity : contains the result of gene variable analysis.
* 02b_FilterData : contains the result of the pre process data used for scBERT.
* 03_scGPTAnalysis : contains the results of the scGPT tool, i.e. cell type predictions.
* 04_scGPTAnalysisResult : contains analysis results of 1 run from the scGPT tool.
* 05_scGPTDifferentEpoch : contains analysis results of multiple runs from the scGPT tool.
* 06_scBERT : contains the results of the scBERT tool, i.e. cell type predictions.

To download and uncompress the data, use the following code:

**On linux:**

```
    cd $WORKING_DIR
    wget https://zenodo.org/records/11658106/files/Human_Thymus_Development_Atlas_output.tar.gz -O Human_Thymus_Development_Atlas_output.tar.gz
    tar zxvf Human_Thymus_Development_Atlas_output.tar.gz
    
    wget https://zenodo.org/records/11658106/files/Human_Thymus_Development_Atlas_Reference.tar.gz -O Human_Thymus_Development_Atlas_Reference.tar.gz
    tar zxvf Human_Thymus_Development_Atlas_Reference.tar.gz
```
 
Once done, you may obtain the following subfolder structure, each of them containing several files.

```
    scTransformers
    ├── cross_tissue_immune_cell
    │   ├── 01_Reference
    │   │   ├── 00_Dataset
    │   │   ├── 01_CellCycle
    │   │   ├── 02_HeatShock
    │   │   └── 04_Model
    │   └── 05_Output
    │       ├── 01_Datapreprocessing
    │       ├── 02a_GlobalHeterogenity
    │       ├── 02b_FilterData
    │       ├── 03_scGPTAnalysis
    │       ├── 04_scGPTAnalysisResult
    │       ├── 05_scGPTDifferentEpoch
    │       ├── 06_scBERT
    │       └── fig
    └── Human_Thymus_Development_Atlas
        ├── 01_Reference
        │   ├── 00_Dataset
        │   ├── 01_CellCycle
        │   ├── 02_HeatShock
        │   └── 04_Model
        └── 05_Output
            ├── 01_Datapreprocessing
            ├── 02a_GlobalHeterogenity
            ├── 02b_FilterData
            ├── 03_scGPTAnalysis
            ├── 04_scGPTAnalysisResult
            ├── 05_scGPTDifferentEpoch
            ├── 06_scBERT
            └── fig

```

---

### Install Singularity and Docker

You need to install Singularity v3.7 on your system to run the complete analysis. Follow the instructions here : https://sylabs.io/guides/3.7/admin-guide/

You also need to install Docker on your system to take advantage of interactive analysis environment with Rstudio and jupyter lab, follow the instructions here : https://docs.docker.com/get-docker/

---

### Download the Singularity images

Singularity images files are stored on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11658106.svg)](https://doi.org/10.5281/zenodo.11658106) (also on Cross tissue immune cell). Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file and untar it:

**On linux:**

```
    cd $WORKING_DIR
    wget https://zenodo.org/records/11658106/files/Human_Thymus_Development_Atlas_container.tar.gz -O Human_Thymus_Development_Atlas_container.tar.gz
    tar zxvf Human_Thymus_Development_Atlas_container.tar.gz
```

These commands will create a sub-folder named **02_Container** in the first dataset folder:

```
    scTransformers
    └── cross_tissue_immune_cell
        └── 02_Container
```

This folder contains the Singularity images for the single-cell RNA-seq analysis. Since the Singularity images are used for the 2 tools, they must be present in all the sample folder in the same **02_Container** subfolder. Instead of copying the image files, we will create symbolic links to spare disk space:

**On linux:**

```
    cd $WORKING_DIR
    ln -s $WORKING_DIR/cross_tissue_immune_cell/02_Container/cb2m_scbert_gpur_27_03_2024.img $WORKING_DIR/Human_Thymus_Development_Atlas/02_Container/cb2m_scbert_gpur_27_03_2024.img
    ln -s $WORKING_DIR/cross_tissue_immune_cell/02_Container/scgpt_test_rf.img $WORKING_DIR/Human_Thymus_Development_Atlas/02_Container/scgpt_test_rf.img
```

---

### Download the Docker images and load them on your system

Docker image tar files are stored on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11658106.svg)](https://doi.org/10.5281/zenodo.11658106) (also on Cross tissue immune cell). Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file, untar it and load the docker images on your system: 

```
    cd $WORKING_DIR
    wget https://zenodo.org/records/11658106/files/Human_Thymus_Development_Atlas_container.tar.gz -O Human_Thymus_Development_Atlas_container.tar.gz
    tar zxvf Human_Thymus_Development_Atlas_container.tar.gz
    docker load -i Human_Thymus_Development_Atlas_container.tar
```

---
---

## Run the analysis

### Run the analysis workflow

The analysis workflow uses the Singularity images and docker pre and post process.

The study contains 2 data. Each data have several steps of analysis for which you will find the R and python script files in the subfolder **03_Script**.

Each step of analysis can generates its own HTML report file and several output files. Some output files of some steps are used by other steps, making a complete workflow of analysis. The output files are stored in each datatset folder in a sub-folder named "05_Output".

### Run the analysis individually using Docker

If you have loaded the docker images (see above), you can use Rstudio or jupyter lab in Docker to run the analysis individually.

To start a docker container, use the following command:

```
docker run -d -p 8787:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  <IMAGE_NAME>
```

where:

* <PASSWORD> is a simple string you will use as password to login into Rstudio
* <IMAGE_NAME> is the Docker image name to run

One started, you can open an internet browser and use the URL https://127.0.0.1:8787.

At the login prompt, enter the name of the user session you are connected with and the password you type in place of <PASSWORD>. You are now in a Rstudio environment and the container is able to connect to the **WORKING_DIR**
of your system. Inside you will find the project files. (Do the same procedure for jupyter lab).

To run the analysis, follow the instruction :

**NOTE** : The tools can use some same script (like variable gene), but they also have their own steps, so please follow carefully :

---
---

## Tools Transformers

---

### scBERT 
![Python 3.8+](https://img.shields.io/badge/python-3.6.8-brightgreen) 

Below you will find detailed instruction of each step. Each step is named after the directory assigned to it. With the objective, and what we should get out of it.

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

**Output:** 
- HTML file by Quarto

---
---

### scGPT 
![Python 3.9+](https://img.shields.io/badge/python-3.6.8-brightgreen) 

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
