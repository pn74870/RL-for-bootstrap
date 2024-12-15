import numpy as np
from scipy.special import hyp2f1
import math
import cmath
from joblib import Parallel, delayed

def kronecker_delta(h, hb, tol=1e-8):
    return 1.0 if abs(h - hb) < tol else 0.0

def g(h, hb, z):
    # Import inside the function if still problematic:
    # from scipy.special import hyp2f1
    delta = kronecker_delta(h, hb)
    denominator = delta + 1.0
    z_conj = np.conjugate(z)

    try:
        hg_h_z = hyp2f1(h, h, 2*h, z)
        hg_hb_z_conj = hyp2f1(hb, hb, 2*hb, z_conj)
        hg_h_z_conj = hyp2f1(h, h, 2*h, z_conj)
        hg_hb_z = hyp2f1(hb, hb, 2*hb, z)
    except Exception:
        return 0.0

    term1 = (z**h)*(z_conj**hb)*hg_h_z*hg_hb_z_conj
    term2 = (z_conj**h)*(z**hb)*hg_h_z_conj*hg_hb_z
    result = (term1 + term2)/denominator
    return result.real

def sample_z(center, width):
    real_part = np.random.normal(loc=center.real, scale=width)
    imag_part = np.random.normal(loc=center.imag, scale=width)
    return complex(real_part, imag_part)

def compute_equation(d_values, z_list, s, Delta_sigma):
    num_equations = len(z_list)
    n = len(d_values)
    A = np.zeros((num_equations, n), dtype=np.float64)
    D = np.zeros(num_equations, dtype=np.float64)

    for eq_idx, z in enumerate(z_list):
        abs_z_minus_1_pow = abs(z - 1)**(2*Delta_sigma)
        abs_z_pow = abs(z)**(2*Delta_sigma)
        D[eq_idx] = - (abs_z_minus_1_pow - abs_z_pow)

        for i in range(n):
            Delta_i = d_values[i]
            s_i = s[i]
            h_i = (Delta_i + s_i)/2
            hb_i = (Delta_i - s_i)/2

            g_z = g(h_i, hb_i, z)
            g_1_minus_z = g(h_i, hb_i, 1 - z)
            coefficient = abs_z_minus_1_pow*g_z - abs_z_pow*g_1_minus_z
            A[eq_idx, i] = coefficient

    return A, D

def compute_C_for_sample(d_values, n_states, s, Delta_sigma, centers, gauss_stds):
    # Ensure gauss_stds is a list if it's scalar:
    if np.isscalar(gauss_stds):
        gauss_stds = [gauss_stds]*4

    z_samples = [sample_z(centers[i % 4], gauss_stds[i % 4]) for i in range(n_states)]
    A, D = compute_equation(d_values, z_samples, s, Delta_sigma)
    try:
        C = np.linalg.solve(A, D)
        return C
    except np.linalg.LinAlgError:
        return None

def parallel_compute_cs_list(d_values, n_states, s, Delta_sigma, centers, gauss_stds, num_samples, n_jobs=-1):
    return Parallel(n_jobs=n_jobs)(
        delayed(compute_C_for_sample)(d_values, n_states, s, Delta_sigma, centers, gauss_stds)
        for _ in range(num_samples)
    )
