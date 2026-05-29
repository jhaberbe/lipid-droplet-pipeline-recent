import sys
import pickle
import tifffile
import itertools
from tqdm import tqdm 
from src.weka import *

folder = sys.argv[1] 

# Load Classifier from weka-type-segmentation
model = pickle.load(open("/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/models/classifier.pickle", "rb"))

# Load Image
img = tifffile.imread(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/slides/plin2/{folder}.tif")

# Normalize
img = normalize_image(img)

# Create logits array
new_img = np.zeros_like(img)

# Iterate over slide 
for i, j in tqdm(list(itertools.product(range((img.shape[0] // 800)), range((img.shape[1] // 800))))):
    try:
        patch = convert_image_array(img[800*i:800*(i+1), 800*j:800*(j+1)])
        logits = model.predict_proba(patch)[:, 1].reshape(800, 800)
        new_img[800*i:800*(i+1), 800*j:800*(j+1)] = logits
    except:
        Warning(f"failed patch at x: {i*800}, y: {j * 800}")

# Saving twice.
pickle.dump(new_img, open(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/segmentation/plin2/{folder}.pickle", "wb"))
tifffile.imwrite(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/segmentation/plin2/{folder}.tif", new_img)