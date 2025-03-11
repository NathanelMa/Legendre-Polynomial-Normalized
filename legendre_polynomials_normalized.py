import numpy as np
import matplotlib.pyplot as plt
from typing import Callable

"""
This script performs a projection of a function onto a set of orthogonal polynomials (Legendre polynomials).
The main functions are as follows:
- `integration_simpson_sampled`: Numerically integrates a function using the Simpson's rule with sampled points.
- `get_projection`: Projects a function onto a set of Legendre polynomials using numerical integration.
- `legendre_polynomial_normalized`: Computes a set of normalized Legendre polynomials over a given interval [a, b].

The Legendre polynomials are used as a basis to approximate a function by projecting it onto this set.
The script also demonstrates the projection of the functions sin(x) and exp(sqrt(abs(x/7))) onto the basis.
"""

def integration_simpson_sampled(a, b, fx, N):
    """
        Numerically integrates a function using the Simpson's rule with sampled points.

        Args:
            a (float): Lower bound of the interval.
            b (float): Upper bound of the interval.
            fx (np.array): Function values at the sampled points.
            N (int): Number of sub intervals to divide the interval [a, b].

        Returns:
            float: The result of the integration.
        """

    h = (b - a) / N
    fx[1:-1:2] *= 4
    fx[2:-1:2] *= 2
    return (h / 3) * np.sum(fx)


def get_projection(BASE: np.array, f: Callable, a: float, b: float):
    """
        Projects a function onto a set of orthogonal polynomials using numerical integration.

        Args:
            BASE (np.array): A 2D array where each row corresponds to a Legendre polynomial.
            f (Callable): The function to be projected.
            a (float): Lower bound of the interval.
            b (float): Upper bound of the interval.

        Returns:
            np.array: The projection of the function onto the polynomials.
        """

    N = BASE.shape[1]
    xx = np.linspace(a, b, N)
    fx = f(xx)
    projection = np.zeros_like(xx)
    for k in range(BASE.shape[0]):
        # sum of: <f,ek> ek
        ek = BASE[k, :]
        projection += integration_simpson_sampled(a, b, fx * ek, N) * ek
        # plt.plot(xx, projection)
        # plt.grid()
        # plt.pause(0.1)

    return projection


def legendre_polynomial_normalized(number_functions, N: int, a, b):
    """
    Generates and normalizes a set of Legendre polynomials over the interval [a, b].

    Args:
        number_functions (int): The number of Legendre polynomials to generate.
        N (int): The number of sampled points to use for numerical integration.
        a (float): The lower bound of the interval.
        b (float): The upper bound of the interval.

    Returns:
        np.array: A 2D array where each row is a normalized Legendre polynomial.
    """
    xx = np.linspace(a, b, N)
    BASE = np.zeros((number_functions, N))
    for n in range(number_functions):
        BASE[n, :] = np.power(xx, n)

    norm_first_poly = integration_simpson_sampled(a, b, np.ones_like(xx), N)
    e1 = np.ones_like(xx) / np.sqrt(norm_first_poly)
    BASE[0, :] = e1

    for n in range(1, number_functions):
        wn = np.power(xx, n)
        p_sum = np.zeros_like(wn)

        for k in range(0, n):
            # <v_n, e_k> e_k
            ek = BASE[k, :]
            p_sum += integration_simpson_sampled(a, b, wn * ek, N) * ek

        wn -= p_sum
        norm = integration_simpson_sampled(a, b, np.power(wn, 2), N)
        BASE[n, :] = wn / np.sqrt(norm)

    return BASE


"""     Legendre Polynomials    """


def main():
    """
    Main function that demonstrates the projection of the functions sin(x) and exp(sqrt(abs(x/7)))
    onto a set of normalized Legendre polynomials and plots the results.

    The projection is computed using the `get_projection` function and the results are displayed alongside
    the original functions.

    The script uses a set of 11 normalized Legendre polynomials over the interval [-pi, pi].
    """
    eps = 1e-10
    number_functions = 11
    a = -np.pi
    b = np.pi
    N = 100
    xx = np.linspace(a, b, N)
    BASE = legendre_polynomial_normalized(number_functions, N, a, b)

    # For plotting all generated polynomial functions (normalized)

    # for i in range(number_functions):
    #   plt.plot(xx, BASE[i])

    # Using the projection over BASE and plot the res
    f = lambda x: np.sin(x)
    g = lambda x: np.exp(np.sqrt(np.abs((1 / 7) * x)))
    projection_f = get_projection(BASE, f, a, b)
    projection_g = get_projection(BASE, g, a, b)

    plt.figure(figsize=(10, 6))
    plt.plot(xx, g(xx), color='darkblue', alpha=0.5, lw=5, linestyle='--',
             label=r'$g(x) = e^{\sqrt{\left|\frac{x}{7}\right|}}$')
    plt.plot(xx, f(xx), color='blue', alpha=0.5, lw=5, linestyle='--', label=f"$f_x=sin(x)$")
    plt.plot(xx, projection_f, color='red', alpha=0.8, lw=2.5, label=f"$f_x projection$")
    plt.plot(xx, projection_g, color='darkred', alpha=0.9, lw=2.5, label=f"$g_x projection$")
    plt.axhline(0, color='black', linewidth=1), plt.axvline(0, color='black', linewidth=1)
    plt.grid(), plt.legend(), plt.legend(fontsize='xx-large'), plt.show()


if __name__ == '__main__':
    main()
