# HST Photometry Pipeline for NGC 1261

A Python pipeline for processing Hubble Space Telescope (HST) WFPC2 images of the globular cluster NGC 1261 to create a Hertzsprung-Russell (HR) diagram.

## Project Description

This project implements a complete photometric pipeline that processes HST observations through cosmic ray removal, star detection, source matching, aperture photometry, and HR diagram creation. The pipeline analyzes images taken with two filters (F336W and F555W) to determine stellar colors and magnitudes, revealing the stellar populations within the globular cluster.

### Features

- **Cosmic Ray Removal**: Median-combines multiple exposures per filter
- **Star Detection**: Uses local maxima detection with 2D Gaussian fitting
- **Quality Filtering**: Applies flux, ellipticity, and size checks to reject spurious detections
- **Cross-Matching**: Matches sources between filters using KD-tree nearest-neighbor search
- **Aperture Photometry**: Circular apertures with local background subtraction
- **HR Diagram**: Creates publication-quality color-magnitude diagram

## Scientific Context

NGC 1261 is a globular cluster containing thousands of stars. By plotting their colors versus magnitudes in an HR diagram, we can identify distinct stellar populations including the main sequence, red giant branch, horizontal branch, and potentially blue stragglers.

The two filters used capture different wavelength ranges:
- **F336W**: Ultraviolet (~336 nm) - sensitive to hot, blue stars
- **F555W**: Visible green (~555 nm) - provides reference magnitude

The color difference (F336W - F555W) serves as a proxy for stellar temperature, while F555W magnitude represents stellar brightness.

## Usage Instructions

Simply run the main script:

```bash
python hst_photometry_pipeline.py
```

The pipeline will automatically:
1. Process all FITS files in `Data/F336W` and `Data/F555W`
2. Generate intermediate and final output files
3. Display the HR diagram


## Output Files

The pipeline generates the following files:

| File | Description |
|------|-------------|
| `F336W_median.fits` | Cosmic-ray-cleaned F336W image |
| `F555W_median.fits` | Cosmic-ray-cleaned F555W image |
| `matched_sources.csv` | Catalog of sources detected in both filters |
| `photometry_catalog.csv` | Full photometry catalog with magnitudes |
| `HR_diagram.png` | Hertzsprung-Russell diagram (300 DPI) |


## Pipeline Details

### Step 1: Cosmic Ray Removal

- **Method**: Median-combining multiple exposures
- **Input**: 3 FITS files per filter
- **Output**: Single median-combined image per filter
- **Rationale**: Cosmic rays are random events; median combination effectively removes them while preserving stellar signals

### Step 2: Star Detection

- **Method**: Local maxima detection + 2D Gaussian PSF fitting
- **Detection threshold**: 0.6σ 
- **Fitting box size**: 9×9 pixels
- **Quality checks**:
  - Positive flux (amplitude > 0)
  - Ellipticity ≤ 0.9 (roundness criterion)
  - Size: 0.3 ≤ σ ≤ 6.0 pixels

### Step 3: Source Matching

- **Method**: KD-tree nearest-neighbor search
- **Matching tolerance**: 2.0 pixels
- **Strategy**: Each source in F336W is matched to the nearest source in F555W if within tolerance
- **Duplicate prevention**: Each F555W source can only match once
- **Output**: ~1,000-1,300 matched sources
- **Match rate**: Typically 85-95%

### Step 4: Aperture Photometry

- **Aperture**: Circular, radius = 3.0 pixels
- **Background**: Local annulus (inner radius = 4.0 pixels, outer radius = 7.0 pixels)
- **Background estimation**: Median of pixels in annulus
- **Flux calculation**: Sum of (aperture pixels - background)
- **Magnitude system**: Instrumental magnitudes
  - Formula: mag = -2.5 × log₁₀(flux)
- **Invalid measurements**: Removed (negative flux, edge sources)

### Step 5: HR Diagram

- **X-axis**: Color (mag_F336W - mag_F555W)
- **Y-axis**: Magnitude (mag_F555W, inverted)
- **Format**: Scatter plot with ~1,000-1,200 stars
- **Convention**: Inverted Y-axis (bright stars at top)
- **Expected features**: Main sequence, red giant branch, horizontal branch

## Assumptions and Limitations

### Assumptions

1. **Image alignment**: F336W and F555W images are pre-aligned to within a few pixels
2. **FITS structure**: Science data is in HDU extension 1 (not PRIMARY HDU)
3. **Stellar PSF**: Approximated as 2D Gaussian (valid for HST WFPC2)
4. **Cosmic ray distribution**: Random and uncorrelated between exposures
5. **Background variation**: Varies smoothly across the field; local measurement is appropriate
6. **Point sources**: All detected objects are assumed to be stars (not extended sources like galaxies)

### Limitations

1. **Single chip**: Uses only WFPC2 PC chip (800×800 pixels), not all 4 chips
2. **Crowded fields**: May miss stars in very dense regions or reject blended sources
3. **Edge effects**: Stars within 5 pixels of image edges are excluded from photometry
4. **Magnitude calibration**: Produces instrumental magnitudes (not calibrated to standard photometric system)
5. **Filter coverage**: Limited to two filters (F336W, F555W); more filters would provide better color information
6. **Photometric accuracy**: Simple aperture photometry; PSF photometry would be more accurate in crowded fields
7. **Detection completeness**: Not complete at faint magnitudes (completeness analysis not performed)

### Known Issues

- **OptimizeWarning**: Some Gaussian fits fail (covariance matrix not estimated) - these sources are automatically rejected
- **Detection asymmetry**: F555W typically detects fewer stars than F336W due to different filter sensitivities
- **Blended sources**: Close stellar pairs may be detected as single elongated sources and rejected by ellipticity check

## Algorithm Parameters

The following parameters can be adjusted in the code if needed:

| Parameter | Default | Location | Description | Effect if Changed |
|-----------|---------|----------|-------------|-------------------|
| `threshold_sigma` | 0.6 (F336W), 0.5 (F555W) | `find_stars()` | Detection threshold (σ above mean) | Lower = more detections, more false positives |
| `box_size` | 9 | `find_stars()` | Size of cutout for Gaussian fitting | Larger = better for extended sources, slower |
| `ellipticity_tolerance` | 0.9 | `find_stars()` | Maximum allowed ellipticity | Higher = accepts more elongated sources |
| `min_sigma` | 0.3 | `find_stars()` | Minimum Gaussian width (pixels) | Lower = accepts smaller/sharper sources |
| `max_sigma` | 6.0 | `find_stars()` | Maximum Gaussian width (pixels) | Higher = accepts larger/extended sources |
| `tolerance` | 2.0 | `match_sources()` | Matching distance tolerance (pixels) | Higher = more matches, more false matches |
| `aperture_radius` | 3.0 | `measure_photometry()` | Photometric aperture radius (pixels) | Larger = more flux, more sky contamination |
| `sky_inner` | 4.0 | `measure_photometry()` | Background annulus inner radius (pixels) | Must be > aperture_radius |
| `sky_outer` | 7.0 | `measure_photometry()` | Background annulus outer radius (pixels) | Larger = more background pixels, smoother estimate |

### Recommended Parameter Ranges

- **Detection threshold**: 0.5-2.0σ (lower for deep detection, higher for cleaner sample)
- **Matching tolerance**: 1.0-3.0 pixels (depends on image alignment quality)
- **Aperture radius**: 2.0-5.0 pixels (should enclose ~90% of stellar flux)
- **Ellipticity tolerance**: 0.5-0.9 (0.5 for very round sources only, 0.9 more permissive)
