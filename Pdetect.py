"""
Pdetect: http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002419 (Sharp & Agrawal)
"""

from mpileLine import mpileLine
from mpSample import mpSample
from scipy.stats.distributions import binom

def BB(i, n, freq):
    """
    """
    pass

def Pdetect(mp_sample):
    """ Calculates the power to detect mutations for a single individual at one site

    Args:
        mp_sample: A mpSample object representing a single sample

    Returns:
        p_detect: float 
    """
    min_num_reads = 5 #Temporary
    f_mut = 0.5

    n_f = mp_sample.get_num_fwd_reads()
    n_r = mp_sample.get_num_rev_reads()

    X_i_j = 0.0
    if n_f + n_r => min_num_reads:
        X_i_j = 1.0

    p_detect = 0.0

    for i in range(1, n_f + 1):
        for j in range(1, n_r + 1):
#            p_detect = p_detect + (X_i_j) * BB(i, n_f, f_mut) * BB(j, n_r, f_mut)
            p_detect = p_dect + (X_i_j) * binom.pmf(i, n_f, f_mut) * binom.pmf(j, n_r, f_mut)

    return p_detect

















