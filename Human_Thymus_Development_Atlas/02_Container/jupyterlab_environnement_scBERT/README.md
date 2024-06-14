# Image for testing Python spatial transcriptomics libraries and interactive data exploration

## Build


### Build from Docker

```
docker build -t <Nom de l'image> [chemin d'accès]
```




## Run

### Run from Docker image

```
docker run -u $(id -u) -e JUPYTER_TOKEN="myPass" -v /mnt/DOSI:/mnt/DOSI -v /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/thyarion/scLLM/:/workspace --name 29_02_2024_test_des_data -p 7777:8888 <nom_image>
```
And connect to port 7777


