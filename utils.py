

import numpy as np
import math

# objective functions
def residual_square(pred_arr, gt_arr):
    return np.sum(np.sqrt(np.square((pred_arr-gt_arr))))

def related_error_rate(pred_arr, gt_arr):
    real_num = pred_arr - gt_arr
    related_errors = real_num / gt_arr
    related_errors_abs = abs(related_errors)

    return np.mean(related_errors_abs)

