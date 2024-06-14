FROM jupyter/minimal-notebook:python-3.8


RUN python -m pip install --upgrade pip

#### Requierement from scBERT
#https://github.com/TencentAILabHealthcare/scBERT/blob/master/requirements.txt
RUN pip install torch
RUN pip install torchvision
RUN pip install transformers==4.6.1
RUN pip install scanpy==1.7.2
RUN pip install scikit-learn==0.24.2
RUN pip install scipy==1.5.4
RUN pip install numpy==1.22
RUN pip install pandas==1.1.5
RUN pip install einops==0.6.0
RUN pip install matplotlib==3.6.0
RUN pip install performer_pytorch

USER root

RUN apt-get install -y git

RUN git clone https://github.com/TencentAILabHealthcare/scBERT /opt/scBERT


ENV NUMBA_CACHE_DIR=/tmp/numba_cache

RUN chmod -R 777 /home/jovyan


# Expose the JupyterLab port
EXPOSE 8888

WORKDIR "/"

# Set the default command to launch JupyterLab
CMD ["jupyter-lab", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root", "--NotebookApp.iopub_data_rate_limit=1.0e10"]
