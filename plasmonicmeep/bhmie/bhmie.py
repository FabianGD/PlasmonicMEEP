import numpy as np


def bhmie(x, refrel, n_angles):
    """
    This file is converted from mie.m, see http://atol.ucsd.edu/scatlib/index.htm
    Bohren and Huffman originally published the code in their book on light scattering

    Literature:
    - https://doi.org/10.1364/AO.43.001951

    Calculation based on Mie scattering theory
    input:
         x      - size parameter = k*radius = 2pi/lambda * radius
                      (lambda is the wavelength in the medium around the scatterers)
         refrel - refraction index (n in complex form for example:  1.5+0.02*i;
         n_angles   - number of angles for S1 and S2 function in range from 0 to pi/2
    output:
           S1, S2 - funtion which correspond to the (complex) phase functions
           Qext   - extinction efficiency
           Qsca   - scattering efficiency
           Qback  - backscatter efficiency
           gsca   - asymmetry parameter
    """

    nmxx = 150000

    s1_1 = np.zeros(n_angles, dtype=np.complex128)
    s1_2 = np.zeros(n_angles, dtype=np.complex128)
    s2_1 = np.zeros(n_angles, dtype=np.complex128)
    s2_2 = np.zeros(n_angles, dtype=np.complex128)

    tau = np.zeros(n_angles)

    if n_angles > 1000:
        print("error: n_angles > max_n_angles=1000 in bhmie")
        return

    # Require NANG>1 in order to calculate scattering intensities
    if n_angles < 2:
        n_angles = 2

    y = x * refrel
    # ymod = abs(m * x)
    ymod = np.abs(y)

    # Series expansion terminated after NSTOP terms
    # Logarithmic derivatives calculated from NMX on down
    # Definition of the break criterion
    criterion_flt = x + 4.0 * x ** 0.3333 + 2.0
    nr_expansion = int(criterion_flt)

    nmx = np.max([criterion_flt, ymod]) + 15.0
    nmx = np.fix(nmx)

    # BTD experiment 91/1/15: add one more term to series and compare results
    #      NMX=AMAX1(XSTOP,YMOD)+16
    # test: compute 7001 wavelengths between .0001 and 1000 micron
    # for a=1.0 micron SiC grain.  When NMX increased by 1, only a single
    # computed number changed (out of 4*7001) and it only changed by 1/8387
    # conclusion: we are indeed retaining enough terms in series!

    if nmx > nmxx:
        print("error: nmx > nmxx=%f for |m|x=%f" % (nmxx, ymod))
        return

    delta_angle = 0.5 * np.pi / (n_angles - 1)
    thetas = np.arange(n_angles, dtype=np.float64) * delta_angle
    amu = np.cos(thetas)

    pi_0 = np.zeros(n_angles)
    pi_1 = np.ones(n_angles)

    # Logarithmic derivative D_n(J) calculated by downward recurrence
    # beginning with initial value (0.,0.) at J=NMX

    max_order = int(nmx) - 1
    logderiv_d = np.zeros(max_order + 1)
    for idx in range(max_order):
        # Assign the order and the index
        order = nmx - idx
        order_idx = max_order - idx

        # Actually compute the log derivate.
        logderiv_d[order_idx - 1] = (order / y) - (
            1.0 / (logderiv_d[order_idx] + order / y)
        )

    # *** Riccati-Bessel functions with real argument X
    #    calculated by upward recurrence

    # In the loop, the following holds:
    # for given N, PSI  = psi_n        CHI  = chi_n
    #              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
    #              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}

    psi_0 = np.cos(x)
    psi_m1 = np.sin(x)
    chi_0 = -np.sin(x)
    chi_m1 = np.cos(x)
    xi_1 = psi_m1 - chi_m1 * 1j
    qsca = 0.0
    gsca = 0.0
    sign = -1

    a_n: float = 0.0
    b_n: float = 0.0
    a_n1: float = 0.0
    b_n1: float = 0.0

    # This loop computes a Taylor series expansion
    for n in range(0, nr_expansion):
        curr_n = n + 1.0
        factor_n = (2.0 * curr_n + 1.0) / (curr_n * (curr_n + 1.0))

        # Calculate psi_n and chi_n
        psi = (2.0 * curr_n - 1.0) * psi_m1 / x - psi_0
        chi = (2.0 * curr_n - 1.0) * chi_m1 / x - chi_0
        xi = psi - chi * 1j

        # *** Store previous values of AN and BN for use
        #     in computation of g=<cos(theta)>
        if n > 0:
            a_n1 = a_n  # type: ignore
            b_n1 = b_n  # type: ignore

        # *** Compute AN and BN:
        a_n = (logderiv_d[n] / refrel + curr_n / x) * psi - psi_m1
        a_n /= (logderiv_d[n] / refrel + curr_n / x) * xi - xi_1

        b_n = (refrel * logderiv_d[n] + curr_n / x) * psi - psi_m1
        b_n /= (refrel * logderiv_d[n] + curr_n / x) * xi - xi_1

        # *** Augment sums for Qsca and g=<cos(theta)>
        qsca += (2.0 * curr_n + 1.0) * (abs(a_n) ** 2 + abs(b_n) ** 2)
        gsca += ((2.0 * curr_n + 1.0) / (curr_n * (curr_n + 1.0))) * (
            np.real(a_n) * np.real(b_n) + np.imag(a_n) * np.imag(b_n)
        )

        if n > 0:
            gsca += ((curr_n - 1.0) * (curr_n + 1.0) / curr_n) * (
                np.real(a_n1) * np.real(a_n)
                + np.imag(a_n1) * np.imag(a_n)
                + np.real(b_n1) * np.real(b_n)
                + np.imag(b_n1) * np.imag(b_n)
            )

        # *** Now calculate scattering intensity pattern
        #     First do angles from 0 to 90
        pi = pi_1.copy()  # 0+pi1 because we want a hard copy of the values
        tau = curr_n * amu * pi - (curr_n + 1.0) * pi_0
        s1_1 += factor_n * (a_n * pi + b_n * tau)
        s2_1 += factor_n * (a_n * tau + b_n * pi)

        # *** Now do angles greater than 90 using PI and TAU from
        #     angles less than 90.
        #     P=1 for N=1,3,...% P=-1 for N=2,4,...
        #     remember that we have to reverse the order of the elements
        #     of the second part of s1 and s2 after the calculation
        sign = -sign
        s1_2 += factor_n * sign * (a_n * pi - b_n * tau)
        s2_2 += factor_n * sign * (b_n * pi - a_n * tau)

        psi_0 = psi_m1
        psi_m1 = psi

        chi_0 = chi_m1
        chi_m1 = chi

        xi_1 = psi_m1 - chi_m1 * 1j

        # *** Compute pi_n for next value of n
        #     For each angle J, compute pi_n+1
        #     from PI = pi_n , PI0 = pi_n-1
        pi_1 = ((2.0 * curr_n + 1.0) * amu * pi - (curr_n + 1.0) * pi_0) / curr_n
        pi_0 = pi.copy()  # 0+pi because we want a hard copy of the values

    # *** Have summed sufficient terms.
    #     Now compute QSCA,QEXT,QBACK,and GSCA

    # we have to reverse the order of the elements of the second part of s1 and s2
    s1 = np.concatenate((s1_1, s1_2[-2::-1]))
    s2 = np.concatenate((s2_1, s2_2[-2::-1]))
    gsca = 2.0 * gsca / qsca
    qsca = (2.0 / x ** 2) * qsca
    qext = (4.0 / x ** 2) * np.real(s1[0])

    # more common definition of the backscattering efficiency,
    # so that the backscattering cross section really
    # has dimension of length squared
    qback = 4 * (np.abs(s1[2 * n_angles - 2]) / x) ** 2
    # qback = ((abs(s1[2*nang-2])/dx)**2 )/pii  #old form

    return s1, s2, qext, qsca, qback, gsca


if __name__ == "__main__":
    bhmie(x=0.1, refrel=1.5 + 0.5j, n_angles=100)
