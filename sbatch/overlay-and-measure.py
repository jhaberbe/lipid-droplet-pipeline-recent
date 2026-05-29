import sys
import pickle
import pathlib
import tifffile
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from tqdm import tqdm

from PIL import Image, ImageDraw
from shapely.geometry import Polygon

oil_red_o_image_mapping = {
    "13-69": "Slide_2_13-69_05-27.tif",
    "05-27": "Slide_2_13-69_05-27.tif",
    "14-02": "Slide_1_14-02_18-20.tif",
    "10-46": "Slide_4_10-46_13-54_Slide_Scan.tif",
    "18-20": "Slide_1_14-02_18-20.tif",
    "13-54": "Slide_4_10-46_13-54_Slide_Scan.tif",
    "04-44": "0029269_AD_33_04-44_15-27.tif",
    "15-27": "0029269_AD_33_04-44_15-27.tif",
    "14-20": "0029281_ND_33_18-75_14-20.tif",
    "18-75": "0029281_ND_33_18-75_14-20.tif",
    "04-06": "0029282_AD_44_99-15_04-06.tif",
    "99-15": "0029282_AD_44_99-15_04-06.tif"
}

def polygons_to_instance_mask(polygons, image_size):
    """
    Convert a list of polygons to an instance segmentation mask.

    Args:
        polygons: List of polygons, each a list of (x, y) tuples.
        image_size: (width, height) of the output mask.

    Returns:
        A numpy array of shape (height, width) with instance labels.
    """
    width, height = image_size
    mask = Image.new("I", (width, height), 0)  # "I" = 32-bit signed integer pixels
    draw = ImageDraw.Draw(mask)

    for idx, poly_coords in tqdm(enumerate(polygons, start=1)):  # Start instance IDs from 1
        # Draw polygon as filled shape
        draw.polygon(poly_coords, fill=idx)

    return np.array(mask)

sys.path.append(directory)
from src.image_alignment import *

if __name__ == "__main__":
    directory = sys.argv[1]
    folder_key = sys.argv[2]

    raw_slide_location = pathlib.Path(f"{directory}/data/raw/slides/")
    segmentation_slide_location = pathlib.Path(f"{directory}/data/processed/segmentation")
    segmentation_alignments_location = pathlib.Path(f"{directory}/data/alignments")

    boundaries = pd.read_csv(list(pathlib.Path(f"{directory}/data/raw/xenium_runs").glob(f"*{folder_key}*"))[0] / "cell_boundaries.csv.gz")

    polygons = [
        [(x, y) for x, y in zip(df["vertex_x"]/.2125, df["vertex_y"]/.2125)]
        for cell, df in boundaries.groupby("cell_id")
    ]

    plin_path = raw_slide_location / "plin2" / f"{folder_key}.tif"
    plin_image = (tifffile.imread(plin_path) >= threshold_by_folder[folder_key])

    oil_red_o_path = segmentation_slide_location / "oil-red-o" / oil_red_o_image_mapping[folder_key]
    oil_red_o_image = tifffile.imread(oil_red_o_path) >= 0.5 

    # Finds the size of the slide, and returns the minimal bounding box plus a large berth for safety.
    cells = pd.read_csv(list(pathlib.Path(f"{directory}/data/raw/xenium_runs").glob(f"*{folder_key}*"))[0] / "cells.csv.gz")
    x_max, y_max = cells[["x_centroid", "y_centroid"]].max().tolist()
    output_shape = (round(y_max / 0.2125), round(x_max / 0.2125))

    print(f"Warping {folder_key} - Oil Red O")
    oil_red_o_image = apply_transformation_matrix(
        oil_red_o_image, 
        segmentation_alignments_location / "oil-red-o" / f"{folder_key}.csv", 
        output_shape
    ) > 0
    instance_segmentation = polygons_to_instance_mask(polygons, oil_red_o_image.shape[::-1])
    oil_red_o_image = instance_segmentation * oil_red_o_image
    tifffile.imwrite(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/warped_segmentation/oil_red_o/{folder_key}.tif", instance_segmentation * (oil_red_o_image).astype(int))

    print(f"Measurement of {folder_key} Oil Red O")
    oil_red_o_measurements = generate_measurements(oil_red_o_image)
    oil_red_o_measurements.to_csv(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/locations/oil-red-o/{folder_key}.csv")