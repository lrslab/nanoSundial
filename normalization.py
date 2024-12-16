import numpy as np
OUTLIER_THRESH = 5.0
MAX_POINTS_FOR_THEIL_SEN = 1000

def normalize_signal_with_lim(raw_signal, lower_lim=-OUTLIER_THRESH, upper_lim=OUTLIER_THRESH):
    # shift = np.median(signal)
    # # note factor to scale MAD to SD
    # scale = np.median(np.absolute(signal - shift)) * 1.4826
    # norm_signal = (signal - shift) / scale
    shift = np.median(raw_signal)
    scale = np.median(np.abs(raw_signal - shift))
    norm_signal = (raw_signal - shift) / scale

    # norm_signal = c_apply_outlier_thresh(
    #         norm_signal, lower_lim, upper_lim)

    norm_signal = np.where(norm_signal > upper_lim, upper_lim, norm_signal)
    norm_signal = np.where(norm_signal < lower_lim, lower_lim, norm_signal)
    return norm_signal,shift,scale

def normalize_signal(raw_signal):
    # shift = np.median(signal)
    # # note factor to scale MAD to SD
    # scale = np.median(np.absolute(signal - shift)) * 1.4826
    shift = np.median(raw_signal)
    scale = np.median(np.abs(raw_signal - shift))
    norm_signal = (raw_signal - shift) / scale

    # norm_signal = c_apply_outlier_thresh(
    #         norm_signal, lower_lim, upper_lim)
    return norm_signal


def calc_kmer_fitted_shift_scale(prev_shift, prev_scale, r_event_means, r_model_means,):
    """Use robust Theil-Sen estimator to compute fitted shift and scale
    parameters based on read sequence

    Args:
        prev_shift (float): previous shift parameter
        prev_scale (float): previous scale parameter
        r_ref_means (`np.array::np.float64`): expected base signal levels
        r_ref_sds (`np.array::np.float64`): expected base signal level sds


    Returns:
        Sequence-fitted scaling parameters

        1) shift parameter (float)
        2) scale parameter (float)
        3) shift correction factor; multiply by ``prev_scale`` and add to ``prev_shift`` to get ``shift`` (float)
        4) scale correction factor; multiply by ``prev_scale`` to get ``scale`` (float)
    """
    def compute_slopes(r_event_means, r_model_means):
        # despite computing each diff twice this vectorized solution is
        # about 10X faster than a list comprehension approach
        delta_event = r_event_means[:, np.newaxis] - r_event_means
        delta_model = r_model_means[:, np.newaxis] - r_model_means
        return delta_model[delta_event > 0] / delta_event[delta_event > 0]

    n_points = r_model_means.shape[0]
    # potentially sample points for long reads (>1kb)
    if r_model_means.shape[0] > MAX_POINTS_FOR_THEIL_SEN:
        n_points = MAX_POINTS_FOR_THEIL_SEN
        samp_ind = np.random.choice(
            r_model_means.shape[0], n_points, replace=False)
        r_model_means = r_model_means[samp_ind]
        r_event_means = r_event_means[samp_ind]
    # compute Theil-Sen slope estimator
    slope = np.median(compute_slopes(r_event_means, r_model_means))
    inter = np.median(r_model_means - (slope * r_event_means))
    if slope == 0:
        raise RuntimeError('Read failed sequence-based signal ' +
                            're-scaling parameter estimation.')
        # convert to shift and scale parameters (e.g. (obs - shift) / scale)
        scale_corr_factor = 1 / slope
        shift_corr_factor = -inter / slope

    # apply shift and scale values fitted from kmer conditional model
    shift = prev_shift + (shift_corr_factor * prev_scale)
    scale = prev_scale * scale_corr_factor
    return shift, scale, shift_corr_factor, scale_corr_factor

# from nanoCEM._c_helper import c_new_means
# def compute_base_means(all_raw_signal, base_starts):
#     """Efficiently compute new base mean values from raw signal and base start
#     positions
#
#     Args:
#         all_raw_signal (`np.array`): raw nanopore signal obervation values
#         base_starts (`np.array::np.int32`): 0-based base start positions within
#             raw signal
#
#     Returns:
#         `np.array::np.float64` containing base mean levels
#     """
#     return c_new_means(all_raw_signal.astype(np.float64), base_starts)
#
# def normal_new(signal,segs):
#     norm_signal,shift,scale = normalize_signal_with_lim
#     (shift, scale, shift_corr_factor,
#      scale_corr_factor) = calc_kmer_fitted_shift_scale(
#         shift, scale,
#         compute_base_means(norm_signal, segs))
#     # re-normalize signal with new fitted parameters
#     norm_signal = (norm_signal - shift_corr_factor) / scale_corr_factor
