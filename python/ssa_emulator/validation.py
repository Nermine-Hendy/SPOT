"""Image-consistency and measured-frame comparison metrics."""
from pathlib import Path
import hashlib
import json
import numpy as np
from .frigate_conditioning import background_residual, load_fits_image

def robust_sigma_local(x):
    x = np.asarray(x, dtype=np.float64)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    sig = 1.4826 * mad
    if not np.isfinite(sig) or sig <= 0:
        sig = np.nanstd(x)
    return float(sig)

def align_background_median(reference, candidate):
    reference = np.asarray(reference, dtype=np.float64)
    candidate = np.asarray(candidate, dtype=np.float64)
    return candidate - np.nanmedian(candidate) + np.nanmedian(reference)

def ssim_metric(reference, candidate):
    r = np.asarray(reference, dtype=np.float64)
    c = np.asarray(candidate, dtype=np.float64)

    from skimage.metrics import structural_similarity

    data_range = max(
        float(np.nanmax(r)),
        float(np.nanmax(c)),
        1.0,
    ) - min(
        float(np.nanmin(r)),
        float(np.nanmin(c)),
        0.0,
    )

    return float(
        structural_similarity(
            r,
            c,
            data_range=max(data_range, 1.0),
        )
    )

def psnr_metric(reference, candidate):
    r = np.asarray(reference, dtype=np.float64)
    c = np.asarray(candidate, dtype=np.float64)
    mse = np.nanmean((r - c) ** 2)

    if mse <= 0:
        return np.inf

    data_range = max(
        float(np.nanmax(r)),
        float(np.nanmax(c)),
        1.0,
    )

    return float(
        20.0 * np.log10(data_range / np.sqrt(mse))
    )

def pearson_metric(reference, candidate):
    r = np.asarray(reference, dtype=np.float64).ravel()
    c = np.asarray(candidate, dtype=np.float64).ravel()

    ok = np.isfinite(r) & np.isfinite(c)

    if np.count_nonzero(ok) < 2:
        return np.nan

    return float(np.corrcoef(r[ok], c[ok])[0, 1])

def direct_pair_metrics(reference, candidate):
    r = np.asarray(reference, dtype=np.float64)
    c = np.asarray(candidate, dtype=np.float64)
    d = c - r

    return {
        "SSIM": ssim_metric(r, c),
        "PSNR_dB": psnr_metric(r, c),
        "RMSE_DN": float(np.sqrt(np.nanmean(d ** 2))),
        "MAE_DN": float(np.nanmean(np.abs(d))),
        "Pearson_r": pearson_metric(r, c),
    }

def target_integrated_flux(image, support_mask):
    image = np.asarray(image, dtype=np.float64)
    support_mask = np.asarray(support_mask, dtype=bool)

    background_pixels = image[~support_mask]
    background_median = float(np.nanmedian(background_pixels))

    residual = image - background_median
    flux = float(
        np.nansum(
            np.maximum(
                residual[support_mask],
                0.0,
            )
        )
    )

    return flux, background_median

def frigate_reference_metrics(reference, candidate):
    # Raw image statistics are reported directly.
    r = np.asarray(reference, dtype=np.float64)
    c_raw = np.asarray(candidate, dtype=np.float64)

    # RMSE/MAE are background-aligned by one scalar median offset only.
    c_aligned = align_background_median(r, c_raw)
    diff_aligned = c_aligned - r
    diff_raw = c_raw - r

    _, _, _, sig_r = background_residual(
        r.astype(np.float32),
        filter_size=51,
    )
    _, _, _, sig_c = background_residual(
        c_raw.astype(np.float32),
        filter_size=51,
    )

    return {
        "median_dn": float(np.nanmedian(c_raw)),
        "p99_dn": float(np.nanpercentile(c_raw, 99)),
        "residual_sigma_dn": float(sig_c),
        "background_aligned_rmse_dn": float(
            np.sqrt(np.nanmean(diff_aligned ** 2))
        ),
        "background_aligned_mae_dn": float(
            np.nanmean(np.abs(diff_aligned))
        ),
        "raw_rmse_dn": float(
            np.sqrt(np.nanmean(diff_raw ** 2))
        ),
        "raw_mae_dn": float(
            np.nanmean(np.abs(diff_raw))
        ),
        "ssim_to_frigate": ssim_metric(r, c_aligned),
        "psnr_to_frigate_db": psnr_metric(r, c_aligned),
        "pearson_r_to_frigate": pearson_metric(r, c_aligned),
        "residual_sigma_ratio_to_frigate": (
            float(sig_c / sig_r) if sig_r > 0 else np.nan
        ),
    }


def load_external_satsim(path, metadata_path, shape, noise_on):
    """Load registered raw SatSim output; reject scaled display PNGs and resizing."""
    if not path or not metadata_path:
        raise ValueError('Actual SatSim mode requires image and provenance JSON paths')
    p = Path(path)
    meta = json.loads(Path(metadata_path).read_text(encoding='utf-8'))
    if str(meta.get('simulator', '')).lower() != 'satsim' or not meta.get('version'):
        raise ValueError('Provenance must identify simulator=SatSim and its version')
    if meta.get('noise_on') is not bool(noise_on):
        raise ValueError('Provenance noise_on does not match the requested comparison')
    if meta.get('downsample') != 4 or meta.get('units') != 'DN':
        raise ValueError('Actual SatSim images must use the validation working grid (downsample=4), units=DN')
    if not meta.get('scene_config'):
        raise ValueError('Record the SatSim scene_config or a path to that configuration')
    if p.suffix.lower() in ('.fits', '.fit', '.fts'):
        a, _ = load_fits_image(p)
    elif p.suffix.lower() == '.npy':
        a = np.load(p, allow_pickle=False)
    elif p.suffix.lower() == '.npz':
        with np.load(p, allow_pickle=False) as data:
            a = data['Image_DN']
    else:
        raise ValueError('Use raw FITS, NPY or NPZ DN arrays, not a display image')
    if a.shape != tuple(shape) or not np.isfinite(a).all():
        raise ValueError('SatSim image must already be finite and registered on the full working grid')
    meta['image_sha256'] = hashlib.sha256(p.read_bytes()).hexdigest()
    return np.asarray(a, dtype=np.float32), meta
