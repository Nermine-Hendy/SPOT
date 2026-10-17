"""Frigate conditioning utilities for ssa_emulator."""

import json
from pathlib import Path
from astropy.wcs import WCS
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, median_filter


def load_fits_image(fits_path):
    """Loads a FITS primary image array and header."""
    from astropy.io import fits

    fits_path = Path(fits_path)
    with fits.open(fits_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        if data is None and len(hdul) > 1:
            data = hdul[1].data
            header = hdul[1].header

    return np.asarray(data, dtype=np.float32), header


def extract_frigate_metadata(fits_path, image, header):
    """Extracts basic metadata from FITS header."""

    class Metadata:
        pass

    meta = Metadata()
    meta.exposure_s = float(
        header.get("EXPTIME", header.get("EXPOSURE", 0.5))
    )
    meta.height, meta.width = image.shape
    return meta


def load_solved_wcs(wcs_path):
    """Loads Astrometry.net WCS header solution (compatible with header-only fits)."""
    from astropy.io import fits

    with fits.open(wcs_path) as hdul:
        for hdu in hdul:
            if "CRVAL1" in hdu.header or "CTYPE1" in hdu.header:
                return WCS(hdu.header)
    raise ValueError(f"No valid WCS found in {wcs_path}")


def extract_wcs_metadata(wcs, image_shape):
    """Extracts field of view dimensions from solved WCS."""

    class WCSMeta:
        pass

    meta = WCSMeta()
    ny, nx = image_shape
    pix_coords = np.array([[0, 0], [nx, 0], [0, ny], [nx, ny]])
    sky_coords = wcs.pixel_to_world(pix_coords[:, 0], pix_coords[:, 1])

    ra_deg = sky_coords.ra.deg
    dec_deg = sky_coords.dec.deg

    meta.fov_x_deg = float(np.max(ra_deg) - np.min(ra_deg))
    meta.fov_y_deg = float(np.max(dec_deg) - np.min(dec_deg))
    return meta


def downsample_image(image, factor):
    """Block-averages an image array by an integer factor."""
    if factor <= 1:
        return image
    ny, nx = image.shape
    new_y = ny // factor
    new_x = nx // factor
    trimmed = image[: new_y * factor, : new_x * factor]
    return trimmed.reshape(new_y, factor, new_x, factor).mean(axis=(1, 3))


def scale_roi(roi_xyxy, factor):
    """Scales full-resolution ROI coordinates to downsampled resolution."""
    return tuple(int(v // factor) for v in roi_xyxy)


def background_residual(image, filter_size=101):
    """Extracts low-frequency background, residual image, and robust sigma."""
    image = np.asarray(image, dtype=np.float32)
    bg = median_filter(image, size=filter_size)
    residual = image - bg

    med = np.nanmedian(residual)
    mad = np.nanmedian(np.abs(residual - med))
    sigma = float(1.4826 * mad)

    norm_residual = (
        (residual - med) / sigma if sigma > 0 else np.zeros_like(residual)
    )
    return bg, residual, norm_residual, sigma


def robust_limits(array, percentiles=(1, 99.8)):
    """Returns robust vmin and vmax for visualization."""
    return tuple(np.percentile(array[np.isfinite(array)], percentiles))


def project_catalogue_to_detector(
    stars,
    wcs,
    image_shape,
    max_mag=13.5,
    max_stars=None,
    downsample=1,
):
    """Projects catalogue stars into image-plane pixel coordinates via WCS."""
    ra = np.asarray(stars.ra_deg)
    dec = np.asarray(stars.dec_deg)
    mag = np.asarray(stars.mag)

    valid_mask = (mag <= max_mag) & np.isfinite(ra) & np.isfinite(dec)
    ra, dec, mag = ra[valid_mask], dec[valid_mask], mag[valid_mask]

    pix_x, pix_y = wcs.world_to_pixel_values(ra, dec)

    full_ny, full_nx = image_shape
    in_frame = (
        (pix_x >= 0)
        & (pix_x < full_nx)
        & (pix_y >= 0)
        & (pix_y < full_ny)
    )

    df = pd.DataFrame(
        {
            "x": pix_x[in_frame] / downsample,
            "y": pix_y[in_frame] / downsample,
            "mag": mag[in_frame],
            "ra_deg": ra[in_frame],
            "dec_deg": dec[in_frame],
        }
    )

    if max_stars and len(df) > max_stars:
        df = df.sort_values("mag").head(max_stars).reset_index(drop=True)

    return df


def render_star_flux(
    shape, projected_stars, zero_mag_flux=1.0, psf_sigma_px=0.0
):
    """Renders unit-flux map of point stars."""
    grid = np.zeros(shape, dtype=np.float32)
    ny, nx = shape

    for _, row in projected_stars.iterrows():
        ix, iy = int(round(row["x"])), int(round(row["y"]))
        if 0 <= ix < nx and 0 <= iy < ny:
            flux = zero_mag_flux * (10.0 ** (-0.4 * row["mag"]))
            grid[iy, ix] += flux

    if psf_sigma_px > 0:
        grid = gaussian_filter(grid, sigma=psf_sigma_px, mode="nearest")

    return grid


def manual_trails_from_roi(
    manual_df, roi_full_xyxy, downsample=1, wcs=None
):
    """Converts manual local ROI trail measurements to global frame coordinates."""
    x0_full, y0_full, _, _ = roi_full_xyxy
    x0_work = int(x0_full // downsample)
    y0_work = int(y0_full // downsample)

    df = manual_df.copy()

    df["x1_full"] = (df["x1_roi"] * downsample) + x0_full
    df["y1_full"] = (df["y1_roi"] * downsample) + y0_full
    df["x2_full"] = (df["x2_roi"] * downsample) + x0_full
    df["y2_full"] = (df["y2_roi"] * downsample) + y0_full

    # Correct global working coordinate alignment
    df["x1_work"] = df["x1_roi"] + x0_work
    df["y1_work"] = df["y1_roi"] + y0_work
    df["x2_work"] = df["x2_roi"] + x0_work
    df["y2_work"] = df["y2_roi"] + y0_work

    return df


def render_trail_flux(
    shape, trail_df, amplitude=1.0, sigma_px=0, num_samples=100
):
    """Renders line streaks directly onto the global image array."""
    grid = np.zeros(shape, dtype=np.float32)
    ny, nx = shape

    for _, row in trail_df.iterrows():
        x1, y1 = row["x1_work"], row["y1_work"]
        x2, y2 = row["x2_work"], row["y2_work"]

        xs = np.linspace(x1, x2, num_samples)
        ys = np.linspace(y1, y2, num_samples)

        for x, y in zip(xs, ys):
            ix, iy = int(round(x)), int(round(y))
            if 0 <= ix < nx and 0 <= iy < ny:
                grid[iy, ix] += amplitude / num_samples

    if sigma_px > 0:
        grid = gaussian_filter(grid, sigma=sigma_px, mode="nearest")

    return grid


def configure_spot_for_detector_maps(
    p, image_shape, exposure_s=0.5, fov_x_deg=1.0, fov_y_deg=1.0,
    deterministic=False, random_seed=12, bg_mean_dn=120.0,
    bg_sigma_dn=2.97, psf_sigma=0.75,
):
    """Scene-matched gain/PSF/noise calibration expressed in the current SPOT schema.

    Input maps contain working-grid DN-equivalent integrated signal, not flux
    rates. The adapter converts to rates before calling production optics.
    """
    from copy import deepcopy
    p = deepcopy(p)
    gain = float(bg_mean_dn / max(bg_sigma_dn**2, 1e-6))
    if not np.isfinite(gain) or gain <= 0 or exposure_s <= 0:
        raise ValueError('Scene-matched calibration requires positive background, gain and exposure')
    p.Camera.Ny, p.Camera.Nx = map(int, image_shape)
    p.Camera.FOV_deg = [float(fov_x_deg), float(fov_y_deg)]
    p.Camera.Shape = 'rectangular'
    p.Time.Exposure_s = p.Time.Cadence_s = float(exposure_s)
    p.Time.NumExposures = 1
    p.Targets = []
    p.Stars.Enabled = p.Background.Enabled = p.Earth.Enabled = p.Atmosphere.Enabled = False
    p.Radiometry.DetectorCalibration = 'legacy_gain'
    p.Detector.RadiometricGain_e_per_unit = gain
    p.Detector.ConversionGain_e_per_DN = gain
    p.Detector.PoissonEnabled = not deterministic
    p.Detector.DarkEnabled = p.Detector.ReadEnabled = p.Detector.BiasEnabled = False
    p.Detector.ReadNoise_e = p.Detector.Bias_DN = 0.0
    p.Detector.ADCEnabled = False
    p.Detector.FullWellEnabled = True
    p.Detector.PRNUMap = p.Detector.BadPixelMask = []
    p.RandomSeed = int(random_seed)
    p.Optics.Enabled = True
    p.Optics.Transmission = 1.0
    p.Optics.ElementTransmission = []
    p.Optics.VignettingEnabled = False
    p.Optics.PSFEnabled = True
    p.Optics.PSFModel = 'gaussian'
    p.Optics.SigmaCenter_pix = float(psf_sigma)
    p.Output.SaveResults = p.Output.Plot = p.Plot.Save = False
    return p


def apply_spot_to_conditioned_maps(p, point_source_flux, diffuse_background_flux, exposure_s=0.5):
    """Frigate scene conditioning through the production optics and detector.

    The supplied methodology blurs point sources, then adds its already estimated diffuse background.
    Preserve that conditioning order. SPOT's production Gaussian uses its
    MATLAB-aligned finite kernel/zero boundary and detector uses MT19937 streams;
    these backend differences are deliberate and recorded, not tuned away.
    Orbit, visibility and target BRDF are bypassed in this conditioned test.
    """
    from copy import deepcopy
    from .imaging import optics, detector
    point = np.asarray(point_source_flux, dtype=np.float64)
    background = np.asarray(diffuse_background_flux, dtype=np.float64)
    shape = (p.Camera.Ny, p.Camera.Nx)
    if point.shape != shape or background.shape != shape:
        raise ValueError('Conditioned maps must match Camera.Ny/Nx')
    if not np.isfinite(point).all() or not np.isfinite(background).all():
        raise ValueError('Conditioned maps must be finite')
    if np.any(point < 0) or np.any(background < 0):
        raise ValueError('Conditioned maps must be nonnegative')
    if not np.isclose(exposure_s, p.Time.Exposure_s):
        raise ValueError('Exposure argument disagrees with configuration')
    active = np.ones(shape, dtype=bool)
    point_signal, stages = optics(p, point / exposure_s, active)
    diffuse_config = deepcopy(p)
    diffuse_config.Optics.PSFEnabled = False
    diffuse_signal, _ = optics(diffuse_config, background / exposure_s, active)
    result = detector(p, point_signal + diffuse_signal, 1, active)
    return {'Image': result.Image_DN.astype(np.float32), 'DetectorStages': result,
            'OpticalStages': stages, 'Configuration': p,
            'Provenance': {'mode': 'Frigate-conditioned production SPOT optics/detector',
                           'independent_full_forward_validation': False}}


def load_validation_catalogue(path, max_mag=13.5):
    """Read the MAT star catalogue required by the validation workflow."""
    from scipy.io import loadmat
    from types import SimpleNamespace
    m = loadmat(path)
    ra, dec, mag = (np.ravel(m[k]).astype(float) for k in ('starsRA', 'starsDEC', 'starsMAG'))
    ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(mag) & (mag <= max_mag)
    return SimpleNamespace(ra_deg=ra[ok] % 360, dec_deg=dec[ok], mag=mag[ok])
