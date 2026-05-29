#!/bin/bash

set -e  # stop if anything fails

# 1. Create venv if it doesn't exist
if [ ! -d ".venv" ]; then
  python3.9 -m venv .venv
fi

# 2. Activate venv
source .venv/bin/activate

# 3. Upgrade pip tooling
pip install --upgrade pip setuptools wheel

# 4. Install dependencies
pip install \
  numpy \
  pandas \
  tqdm \
  torch \
  tifffile \
  pyometiff \
  transformers \
  pillow \
  shapely \
  scanpy \
  anndata \
  scipy \
  matplotlib \
  torchvision \
  scikit-image \
  markov_random_field \ 
  --only-binary=:all:

echo "✅ Environment setup complete"