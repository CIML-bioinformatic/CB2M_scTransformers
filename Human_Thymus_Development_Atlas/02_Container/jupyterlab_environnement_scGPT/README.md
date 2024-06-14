# Image for testing Python spatial transcriptomics libraries and interactive data exploration

## Build


### Build from Docker

```
docker build -t <Nom de l'image> [chemin d'accès]
```


## Run

### Run from Docker image

```
docker run -u $(id -u) -e JUPYTER_TOKEN="myPass" -v /mnt/DOSI:/mnt/DOSI -v /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/thyarion/scLLM/:/workspace --name 29_02_2024_test_des_data -p 9999:8888 <nom_image>
```
And connect to port 9999

### Singularity
```
singularity run --nv --writable-tmpfs --env JUPYTER_TOKEN="myPass" --env JUPYTERLAB_WORKSPACES_DIR="/" -B /mnt/DOSI:/mnt/DOSI -B /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE:/home/thyarion/workspace /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/thyarion/ciml-cb2mproj/Project/scLLM/scGPT/02_Container/jupyterlab_environnement_scGPT_TESTrf/scgpt_test_rf.sif
```
