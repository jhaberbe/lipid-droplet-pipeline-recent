# Steps 
QOL:
- Make a single script to run the whole thing.
- Parallelize it further (IHC detection.)


## [X] Setup Directories

TODO: Make this part of the overall script
```bash
export PROJECT_HOME="/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline"
```

### Xenium
```bash
cd /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/xenium
sh /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/xenium/setup-files.sh
```

### PLIN2
```bash
cd /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/slides/plin2
sh /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/slides/plin2/setup-files.sh
```

### Oil Red O
```bash
cd /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/slides/oil-red-o
sh /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/slides/oil-red-o/setup-files.sh
```

### Amyloid
```bash
cd /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/slides/oil-red-o
```

## [X] Run PLIN2 Random Forest Classifier
[partial] Tested 

```bash
cd $PROJECT_HOME
sbatch sbatch/random-forest-image-classfier.sh
```

Which will populate the "processed/segmentation/plin2" folder.

## [X] Run Oil-Red-O Classifier
[ ] Tested 

```bash
cd $PROJECT_HOME
sbatch sbatch/segformer-image-classifier.sh
```

Which will populate the "processed/segmentation/oil-red-o" folder.

## [X] Run Double Positive Detection
[X] Tested 

```bash
cd $PROJECT_HOME
sbatch sbatch/overlay-and-measure.sh
```

## [X] Map Pathology Data to adata objects.
[X] Tested 

```bash
cd $PROJECT_HOME
sbatch sbatch/align-pathology-data.sh
```

## Standard Single Cell Pipeline
[ ] Tested
```bash
cd $PROJECT_HOME
sbatch sbatch/single-cell-pipeline.sh
```