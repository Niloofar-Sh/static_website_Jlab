#!/bin/bash

# Install Voila
pip install voila


# Create a custom kernel that uses Voila
python -m ipykernel install --user --name=voila --display-name="Python (Voila)"

# Configure JupyterLab to launch Voila when a notebook is opened
mkdir -p ~/.jupyter/lab/user-settings/@jupyterlab/notebook-extension/
echo '{"voila":{"enabled":true}}' > ~/.jupyter/lab/user-settings/@jupyterlab/notebook-extension/tracker.jupyterlab-settings
