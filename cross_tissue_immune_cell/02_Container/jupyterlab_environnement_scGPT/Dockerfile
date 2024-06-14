FROM jupyter/minimal-notebook:python-3.9

USER root 

RUN python -m pip install --upgrade pip

#### Requierement from scGPT


# Update the package list
RUN apt-get update -y

# Install git
RUN apt-get install -y git

# Install r-base and tzdata
RUN apt-get install -y r-base tzdata

#Install nvcc
RUN conda update -y -c defaults conda
RUN conda install -y --solver=classic conda-libmamba-solver
RUN conda install -y --solver=libmamba -c nvidia cuda


# Install Python packages using pip
RUN pip install torch
RUN pip install packaging
RUN pip install scgpt "flash-attn<1.0.5"
RUN pip install markupsafe==2.0.1
RUN pip install wandb


ENV NUMBA_CACHE_DIR=/tmp/numba_cache

RUN chmod -R 777 /home/jovyan


# Expose the JupyterLab port
EXPOSE 8888

WORKDIR "/"

# Set the default command to launch JupyterLab
CMD ["jupyter-lab", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root", "--NotebookApp.iopub_data_rate_limit=1.0e10"]
