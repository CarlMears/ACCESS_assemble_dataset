#  From O:\amsr_l2\v05\resampling\show_antenna_patterns.f
#
#  These are Frank's comments on the fortran code.
#
#  The values of these AMSR-E parameters come from
#  E:\AMSR\Resampling_weights\PM1_AMSR_MOD.for
#
#  The values of the Midori2 parameters come from
#  E:\AMSR\Resampling_weights\ADEOSII_AMSR_MOD.for
#
#  For these antenna parameters, Peter found one set for 89 GHZ, which was the
#  averaged of 89A and 89B. This value was put in the IFREQ=8 slot The IFREQ=6
#  and 7 slots were set to 0 For Midori2, the IFREQ=6 and 7 slots contain the
#  50.3 and 52.8 patterns respectively.
#
#  Here, I use slot 7 for 89A and slot 8 for 89B, so they have the same value
#  When I do Midori-2, I will need to decide what to do about the 50 GHz
#  channels (maybe put in average value into slot 6)

#  Comments below by CM, Summer 2021
#
#  converted to python so I can see what's going on
#
#  The half-power points did not come close what is published for AMSR2
#  Partly because of the larger antenna for AMSR2, but also partly because
#  even the AMSRE numbers did not match.
#
#  To deal with this, I introduced angle_scale_fact that scale the input delta
#  so that the beam widths match.
#  angle_scale_fact = np.array([0.8721,0.8356,0.8591,0.9786,0.8185])
#
#  The widths now match what is in the JAXA documentation
#                6,7  11  19   24   37
#  jaxa_width = [1.8,1.2,0.65,0.75,0.35]
#  89 has not been addressed yet....

import numpy as np


def AMSR2_antenna_gain(delta, band):
    #       0         1          2         3         4        5       6       7
    #    7GHZ      11GHZ      19GHZ     24GHZ     37GHZ
    # sin_anglim = [0.105, 0.105, 0.0524, 0.0524, 0.0262, 0.000, 0.0262, 0.0262]
    ant_approx_a = np.array(
        [4.343e-6, 2.096e-6, 1.890e-6, 1.623e-6, 7.248e-7, 0.000, 2.070e-6, 2.070e-6],
        dtype=np.float64,
    )
    ant_approx_b = np.array(
        [6.892e-4, 4.059e-4, 3.727e-4, 7.251e-4, 3.051e-4, 0.000, 2.381e-4, 2.381e-4],
        dtype=np.float64,
    )
    ant_approx_c = np.array(
        [0.503, 0.662, 1.391, 1.804, 1.964, 0.000, 4.593, 4.593], dtype=np.float64
    )
    ant_approx_d = np.array(
        [0.651, 1.345, 4.844, 4.721, 15.18, 0.000, 79.785, 79.785], dtype=np.float64
    )

    angle_scale_fact = np.array([0.8721, 0.8356, 0.8591, 0.9786, 0.8185])
    coeff_a = ant_approx_a[band]
    coeff_b = ant_approx_b[band]
    coeff_c = ant_approx_c[band]
    coeff_d = ant_approx_d[band]

    delta = delta / angle_scale_fact[band]

    gain = (
        coeff_a + coeff_b * np.exp(-coeff_c * delta) + np.exp(-coeff_d * delta**2.)
    )

    return gain


def AMSRE_antenna_gain(delta, band):
    #  0         1          2         3         4        5       6     7
    #  7GHZ      11GHZ      19GHZ     24GHZ     37GHZ
    # sin_anglim = [0.105, 0.105, 0.0524, 0.0524, 0.0262, 0.000, 0.0262, 0.0262]
    ant_approx_a = np.array(
        [4.343e-6, 2.096e-6, 1.890e-6, 1.623e-6, 7.248e-7, 0.000, 2.070e-6, 2.070e-6],
        dtype=np.float64,
    )
    ant_approx_b = np.array(
        [6.892e-4, 4.059e-4, 3.727e-4, 7.251e-4, 3.051e-4, 0.000, 2.381e-4, 2.381e-4],
        dtype=np.float64,
    )
    ant_approx_c = np.array(
        [0.503, 0.662, 1.391, 1.804, 1.964, 0.000, 4.593, 4.593], dtype=np.float64
    )
    ant_approx_d = np.array(
        [0.651, 1.345, 4.844, 4.721, 15.18, 0.000, 79.785, 79.785], dtype=np.float64
    )

    # angle_scale_fact = np.array([0.8721, 0.8356, 0.8591, 0.9786, 0.8185])
    coeff_a = ant_approx_a[band]
    coeff_b = ant_approx_b[band]
    coeff_c = ant_approx_c[band]
    coeff_d = ant_approx_d[band]

    # delta = delta/angle_scale_fact[band]

    gain = (
        coeff_a + coeff_b * np.exp(-coeff_c * delta) + np.exp(-coeff_d * delta * delta)
    )

    return gain


def target_gain(delta_km, diameter_in_km=30.0):
    # using the shape from the 11 GHz footprint from above

    coeff_a = 2.096e-6
    coeff_b = 4.059e-4
    coeff_c = 0.662
    coeff_d = 1.345

    scale_factor = 0.047933 * 30.0 / diameter_in_km
    delta = delta_km * scale_factor
    gain = (
        coeff_a + coeff_b * np.exp(-coeff_c * delta) + np.exp(-coeff_d * delta * delta)
    )

    return gain
