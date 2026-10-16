import numpy as np


def validate_recording(full_t, full_phi):
    """Return real, finite 1-D data suitable for uniformly sampled FFT analysis.

    Sampling intervals may differ from their mean by rtol=1e-4, atol=1e-12 s.
    Larger deviations require resampling before reconstruction.
    """
    try:
        full_t = np.asarray(full_t)
        full_phi = np.asarray(full_phi)
    except (TypeError, ValueError) as exc:
        raise ValueError("Recording timestamps and signal must be numeric arrays.") from exc
    if np.iscomplexobj(full_t) or np.iscomplexobj(full_phi):
        raise ValueError("Recording timestamps and signal must be real-valued.")
    try:
        full_t = full_t.astype(float, copy=False)
        full_phi = full_phi.astype(float, copy=False)
    except (TypeError, ValueError) as exc:
        raise ValueError("Recording timestamps and signal must be numeric arrays.") from exc

    if full_t.ndim != 1 or full_phi.ndim != 1 or full_t.shape != full_phi.shape:
        raise ValueError("Recording timestamps and signal must be matching 1-D arrays.")
    if full_t.size < 2:
        raise ValueError("A recording must contain at least two samples.")
    if not np.all(np.isfinite(full_t)) or not np.all(np.isfinite(full_phi)):
        raise ValueError("Recording timestamps and signal must contain only finite values.")

    with np.errstate(over="ignore", invalid="ignore"):
        intervals = np.diff(full_t)
        dt = np.mean(intervals)
    if not np.all(np.isfinite(intervals)) or not np.isfinite(dt):
        raise ValueError("Recording sampling intervals must be finite.")
    if not np.all(intervals > 0):
        raise ValueError("Recording timestamps must be strictly increasing.")
    if not np.allclose(intervals, dt, rtol=1e-4, atol=1e-12):
        raise ValueError(
            "FFT reconstruction requires approximately uniform sampling "
            "(rtol=1e-4, atol=1e-12 seconds); resample this recording first."
        )
    return full_t, full_phi


def combine_data(full_t, full_phi, phi_an):
    """Prepend constant initial-value padding without trimming the recording.

    The FFT length is twice the smallest power of two >= N (minimum 2048).
    ``phi_an`` is retained for call compatibility; padding does not use a fit.
    """
    full_t, full_phi = validate_recording(full_t, full_phi)
    dt = np.mean(np.diff(full_t))
    N = len(full_phi)
    base_length = max(2048, 1 << (N - 1).bit_length())
    N_added = 2 * base_length - N
    with np.errstate(over="ignore", invalid="ignore"):
        t_neg = full_t[0] + np.arange(-N_added, 0, dtype=np.int64) * dt
    phi_neg = np.full(N_added, full_phi[0])

    fourier_t = np.concatenate([t_neg, full_t])
    fourier_phi = np.concatenate([phi_neg, full_phi])
    if not np.all(np.isfinite(fourier_t)) or not np.all(np.diff(fourier_t) > 0):
        raise ValueError("Cannot represent finite, strictly increasing padded timestamps.")
    return fourier_t, fourier_phi, N_added


def remove_noise(full_t, full_phi, threshold):
    
    dy = np.abs(np.diff(full_phi) / np.diff(full_t))  # derivative
    start_index = np.argmax(dy > threshold)

    start_mean = np.mean(full_phi[:start_index])
    if (start_mean < 0):
        full_phi[start_index:] += start_mean
    else:
        full_phi[start_index:]  -= start_mean
    
    for i in range(start_index):
        full_phi[i] = 0

    return full_phi
