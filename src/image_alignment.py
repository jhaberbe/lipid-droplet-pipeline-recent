import pandas as pd
import skimage.measure
import skimage.transform

def apply_transformation_matrix(image, path_to_matrix, output_shape):
    similarity_transformation_matrix = skimage.transform.SimilarityTransform(pd.read_csv(path_to_matrix, header=None))
    image = skimage.transform.warp(image, similarity_transformation_matrix.inverse, output_shape=output_shape)
    return image

def generate_measurements(image):
    region_props = skimage.measure.regionprops(image)

    result = pd.DataFrame(dict(
        x_centroid = [x.centroid[1] for x in region_props],
        y_centroid = [x.centroid[0] for x in region_props],
        area = [x.area for x in region_props]
    ))

    return result