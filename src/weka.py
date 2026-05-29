import pathlib
import PIL.Image
import itertools
import numpy as np
from tqdm import tqdm
from scipy import ndimage
from scipy.ndimage import gaussian_filter

def normalize_image(image):
    image = image.astype(float).copy()
    image -= np.mean(image)
    image /= np.std(image)
    return image

def image_to_patches(image, patch_size=19):
    """
    Converts an image into rows of nxn patches for each pixel.

    Parameters:
        image (numpy.ndarray): Input image as a 2D array.
        patch_size (int): Size of the patch (nxn).

    Returns:
        numpy.ndarray: Array of patches, where each row corresponds to a patch.
    """
    pad_size = patch_size // 2
    padded_image = np.pad(image, pad_size, mode='constant', constant_values=0)
    patches = []

    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            patch = padded_image[i:i + patch_size, j:j + patch_size]
            patches.append(patch.flatten())

    return np.array(patches)

def compute_hessian_matrix(image):
    dxx = ndimage.gaussian_filter(image, sigma=1, order=(2, 0)).reshape(-1, 1)
    dyy = ndimage.gaussian_filter(image, sigma=1, order=(0, 2)).reshape(-1, 1)
    dxy = ndimage.gaussian_filter(image, sigma=1, order=(1, 1)).reshape(-1, 1)

    return np.concatenate([dxx, dxy, dxy, dyy], axis=1)

def compute_hessian_statistics(row):
    """
    Compute various statistics from the Hessian matrix components.

    Parameters:
        a, b, c, d (float): Components of the Hessian matrix.

    Returns:
        dict: A dictionary containing the computed statistics.
    """
    a, b, c, d = row[0], row[1], row[2], row[3]
    t = 13 / 4  # Gamma normalization factor

    # Compute statistics
    module = np.sqrt(a**2 + b * c + d**2)
    trace = a + d
    determinant = a * d - c * b
    first_eigenvalue = (a + d) / 2 + np.sqrt(4 * b**2 + (a - d)**2) / 2
    second_eigenvalue = (a + d) / 2 - np.sqrt(4 * b**2 + (a - d)**2) / 2
    orientation = 0.5 * np.arccos((4 * b**2 + (a - d)**2) / (4 * b**2 + (a - d)**2 + 1e-10))
    # gamma_normalized_square_eigenvalue_diff = t**4 * (a - d)**2 / ((a - d)**2 + 4 * b**2 + 1e-10) # Currently not implemented
    # square_gamma_normalized_eigenvalue_diff = t**2 * ((a - d)**2 + 4 * b**2)                      # Currently not implemented

    return np.array([
        module,
        trace,
        determinant,
        first_eigenvalue,
        second_eigenvalue,
        orientation,
        # gamma_normalized_square_eigenvalue_diff, # Currently not implemented
        # square_gamma_normalized_eigenvalue_diff, # Currently not implemented
    ])

    # NOTE: Consider implementing as a dictionary
    # return {
    #     "module": module,
    #     "trace": trace,
    #     "determinant": determinant,
    #     "first_eigenvalue": first_eigenvalue,
    #     "second_eigenvalue": second_eigenvalue,
    #     "orientation": orientation,
    #     "gamma_normalized_square_eigenvalue_diff": gamma_normalized_square_eigenvalue_diff,
    #     "square_gamma_normalized_eigenvalue_diff": square_gamma_normalized_eigenvalue_diff,
    # }

def convert_image_array(image):
    raw = [image_to_patches(np.asarray(image), patch_size=19)]
    gaussian = [image_to_patches(gaussian_filter(image, sigma=2**i)) for i in range(3, 4)]
    raw_hessian = [compute_hessian_matrix(image)]
    gaussian_hessian = [compute_hessian_matrix(gaussian_filter(image, sigma=2**i)) for i in range(3, 4)]
    raw_hessian_statistics = [np.vectorize(compute_hessian_statistics, signature='(n)->(8)')(x) for x in raw_hessian]
    gaussian_hessian_statistics = [np.vectorize(compute_hessian_statistics, signature='(n)->(8)')(x) for x in gaussian_hessian]
    sobel = [image_to_patches(ndimage.sobel(np.asarray(image)))]

    return np.concatenate(raw + gaussian + raw_hessian + gaussian_hessian + raw_hessian_statistics + gaussian_hessian_statistics + sobel, axis=1)