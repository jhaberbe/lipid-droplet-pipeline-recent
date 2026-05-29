import sys
import pathlib
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing
from torch import nn
import torch.quantization
import tifffile
from pyometiff import OMETIFFReader
import markov_random_field.markov_random_field
from transformers import pipeline
from transformers import AutoImageProcessor, SegformerForSemanticSegmentation

def read_data_chunk(i: int, j: int):
    chunk = img_array[i*800:(i+1)*800, j*800:(j+1)*800]
    inputs = processor(chunk, return_tensors="pt")
    outputs = model(**inputs)
    logits = outputs.logits  # produces shape (batch_size, num_labels, height/4, width/4)

    # rescale logits to original image size (800, 800)
    upsampled_logits = nn.functional.interpolate(
        logits,
        size=(800, 800),
        mode='bilinear',
        align_corners=False
    )

    # TODO: Change to output the logits.
    # apply argmax on the class dimension
    pred_seg = upsampled_logits.argmax(dim=1)[0]

    # update the patch
    result_img_array[i*800:(i+1)*800, j*800:(j+1)*800] = pred_seg.detach().numpy()


slide_number = sys.argv[1]
if __name__ == "__main__":
    # IHC Slides 
    ihc_slide_paths = list(pathlib.Path("/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet/data/raw/ihc-oil-red-o/").glob("*.ome.tif"))
    slide = ihc_slide_paths[slide_number]

    # Pipeline setup for segformer
    pipe = pipeline("image-segmentation", model="jhaberbe/segformer-b0-finetuned-lipid-droplets-v2")
    model = SegformerForSemanticSegmentation.from_pretrained("jhaberbe/segformer-b0-finetuned-lipid-droplets-v2")
    processor = AutoImageProcessor.from_pretrained("jhaberbe/segformer-b0-finetuned-lipid-droplets-v2")

    # Read the files
    reader = OMETIFFReader(fpath=slide)
    img_array, metadata, xml_metadata = reader.read()

            # Initialize an empty array to store the results
    result_img_array = np.empty_like(img_array[:, :, 0])

    # TODO: parallelize this using multiprocessing or Dask
    # for each of the chunks, read and process the data
    num_chunks_x = img_array.shape[0] // 800
    num_chunks_y = img_array.shape[1] // 800
    for i, j in tqdm(list(itertools.product(range(num_chunks_x), range(num_chunks_y)))):
        read_data_chunk(i, j)

    tifffile.imwrite(f'/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/segmentation/oil-red-o/{slide.name.split(".")[0]}.tif', result_img_array)
