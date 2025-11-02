import glob
import numpy as np
from astropy.io import fits
from scipy.ndimage import maximum_filter
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================================
# STEP 1: COSMIC RAY REMOVAL
# ============================================================================

def median_combine(folder: str, output_path: str) -> np.ndarray:
    """ Median-combine all FITS files in a folder to remove cosmic rays."""
    files = sorted(glob.glob(f"{folder}/*.fits"))
    
    images = []
    for f in files:
        with fits.open(f) as hdul:
            images.append(hdul[1].data.astype(np.float32))
    
    stack = np.dstack(images)
    median_img = np.median(stack, axis=2)
    
    fits.PrimaryHDU(median_img).writeto(output_path, overwrite=True)
    return median_img


# ============================================================================
# STEP 2: STAR FINDING
# ============================================================================

def gaussian_2d(xy, amp, xo, yo, sx, sy, offset):
    """ 2D Gaussian model for star fitting."""
    x, y = xy
    g = offset + amp * np.exp(-(((x - xo)**2) / (2 * sx**2) + ((y - yo)**2) / (2 * sy**2)))
    return g.ravel()


def find_stars(image, threshold_sigma=1, box_size=9, edge_margin=40):
    """ Find stars in image using local maxima detection and Gaussian fitting. 
    Returns array of (x_center, y_center) for detected stars.
    """
    # Find local maxima above threshold
    mean, std = np.mean(image), np.std(image)
    local_max = maximum_filter(image, size=3) == image
    detected = (image > mean + threshold_sigma * std) & local_max
    y_peaks, x_peaks = np.where(detected)
    
    stars = []
    half = box_size // 2
    h, w = image.shape
    
    for y0, x0 in zip(y_peaks, x_peaks):
        # Skip sources too close to edges to avoid count noise as stars
        if x0 < edge_margin or y0 < edge_margin:  # Only left (x) and bottom (y)
            continue
        
        # Extract cutout around peak
        y_min, y_max = y0 - half, y0 + half + 1
        x_min, x_max = x0 - half, x0 + half + 1
        
        # Skip edge sources
        if y_min < 0 or x_min < 0 or y_max > h or x_max > w:
            continue
        
        cutout = image[y_min:y_max, x_min:x_max]
        y, x = np.mgrid[0:cutout.shape[0], 0:cutout.shape[1]]
        
        # Fit Gaussian and apply quality checks
        try:
            p0 = [cutout.max() - cutout.min(), box_size/2, box_size/2, 1.5, 1.5, np.median(cutout)]
            popt, _ = curve_fit(gaussian_2d, (x, y), cutout.ravel(), p0=p0)
            amp, xo, yo, sx, sy, _ = popt
            
            # Quality checks: positive flux, roundness, reasonable size
            if amp > 0 and (1 - min(sx, sy) / max(sx, sy)) <= 0.9 and 0.3 <= (sx + sy) / 2 <= 6.0:
                # Convert local to global coordinates
                stars.append([x_min + xo, y_min + yo])
        except:
            continue
    
    return np.array(stars)

def visualize_detections(image, stars, title="Detected Stars"):
    """Show detected stars overlaid on image."""
    vmin = np.percentile(image, 5)
    vmax = np.percentile(image, 99)
    
    plt.figure(figsize=(10, 10))
    plt.imshow(image, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
    plt.scatter(stars[:, 0], stars[:, 1], s=30, facecolors='none', 
                edgecolors='red', linewidths=0.5, alpha=0.8)
    plt.title(f"{title} ({len(stars)} sources)")
    plt.tight_layout()
    plt.savefig(f"{title.replace(' ', '_')}.png", dpi=150)
    plt.show()

# ============================================================================
# STEP 3: SOURCE MATCHING
# ============================================================================

def match_sources(sources1, sources2, tolerance=2.5):
    """Match sources between two catalogs based on proximity."""
    tree = cKDTree(sources2)
    
    matched = []
    idx1 = []
    idx2 = []
    used = set()
    
    for i, pos1 in enumerate(sources1):
        dist, j = tree.query(pos1)
        
        if dist < tolerance and j not in used:
            matched.append((pos1 + sources2[j]) / 2)
            idx1.append(i)
            idx2.append(j)
            used.add(j)
    
    return np.array(matched), np.array(idx1), np.array(idx2)


# ============================================================================
# STEP 4: APERTURE PHOTOMETRY
# ============================================================================

def measure_photometry(image_f336w, image_f555w, matched_coords, aperture_radius=3.0):
    """Measure photometry for all matched sources in both filters."""
    h, w = image_f336w.shape
    yy, xx = np.ogrid[:h, :w]
    
    magnitudes_f336w = []
    magnitudes_f555w = []
    
    for x, y in matched_coords:
        distance = np.sqrt((xx - x)**2 + (yy - y)**2)
        
        # Aperture and background masks
        aperture = distance <= aperture_radius
        sky = (distance >= 4.0) & (distance <= 7.0)
        
        # Measure flux in both filters
        for image, mag_list in [(image_f336w, magnitudes_f336w), (image_f555w, magnitudes_f555w)]:
            if aperture.any() and sky.any():
                background = np.median(image[sky])
                flux = np.sum(image[aperture] - background)
                mag = -2.5 * np.log10(flux) if flux > 0 else np.nan
            else:
                mag = np.nan
            mag_list.append(mag)
    
    # Create catalog
    catalog = pd.DataFrame({
        'Object_ID': np.arange(1, len(matched_coords) + 1),
        'x_center': matched_coords[:, 0],
        'y_center': matched_coords[:, 1],
        'aperture_radius': aperture_radius,
        'mag_F336W': magnitudes_f336w,
        'mag_F555W': magnitudes_f555w
    })
    
    # Remove invalid magnitudes
    catalog = catalog.dropna()
    return catalog


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    # STEP 1: Cosmic Ray Removal
    median_f336w = median_combine("Data/F336W", "F336W_median.fits")
    median_f555w = median_combine("Data/F555W", "F555W_median.fits")
    
    # STEP 2: Star Finding
    stars_f336w = find_stars(median_f336w, threshold_sigma=0.5, edge_margin=45) #skip 45 pixels from the selected limits. Stimated through plotting
    stars_f555w = find_stars(median_f555w, threshold_sigma=0.5, edge_margin=40)
    
    visualize_detections(median_f336w, stars_f336w, "F336W Detected Stars")
    visualize_detections(median_f555w, stars_f555w, "F555W Detected Stars")

    print(f"F336W: {len(stars_f336w)} stars detected")
    print(f"F555W: {len(stars_f555w)} stars detected")
    print()
    
    # STEP 3: Source Matching
    matched_coords, idx_f336w, idx_f555w = match_sources(stars_f336w, stars_f555w, tolerance=2.0)
    
    catalog = pd.DataFrame({
        'Object_ID': np.arange(1, len(matched_coords) + 1),
        'x_center': matched_coords[:, 0],
        'y_center': matched_coords[:, 1]
    })
    catalog.to_csv('matched_sources.csv', index=False)
    
    print(f"Matched stars: {len(matched_coords)}")
    print(f"Match rate: {len(matched_coords)/min(len(stars_f336w), len(stars_f555w))*100:.1f}%")
    print(f"Catalog saved: matched_sources.csv")
    print()
    
    # STEP 4: Aperture Photometry
    photometry_catalog = measure_photometry(median_f336w, median_f555w, matched_coords, aperture_radius=3.0)
    photometry_catalog.to_csv('photometry_catalog.csv', index=False)
    
    print(f"Photometry complete: {len(photometry_catalog)} stars with valid magnitudes")
    print(f"Catalog saved: photometry_catalog.csv")
    print()
    
    # STEP 5: Hertzsprung-Russell Diagram
    color = photometry_catalog['mag_F336W'] - photometry_catalog['mag_F555W']
    magnitude = photometry_catalog['mag_F336W']
    
    plt.figure(figsize=(8, 10))
    plt.scatter(color, magnitude, s=1, c='black', alpha=0.5)
    plt.xlabel('$m_{F336W} - m_{F555W}$', fontsize=14)
    plt.ylabel('$m_{F336W}$', fontsize=14)
    plt.title('NGC 1261 - Hertzsprung-Russell Diagram', fontsize=16)
    plt.gca().invert_yaxis()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('HR_diagram.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"HR diagram created with {len(color)} stars")
    print(f"Saved as: HR_diagram.png")