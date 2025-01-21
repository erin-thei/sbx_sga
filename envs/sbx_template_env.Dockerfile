FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_sga_env

COPY envs/sbx_sga_env.yml ./

# Install environment
RUN conda env create --file sbx_sga_env.yml --name sbx_sga

ENV PATH="/opt/conda/envs/sbx_sga/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_sga", "/bin/bash", "-c"]

# Run
CMD "bash"