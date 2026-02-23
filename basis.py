"""
Basis construction for minimal free resolutions.

Provides functions to construct basis elements for the minimal free
resolution of the binomial edge ideal of a complete graph K_n.
"""

from sage.all import Combinations
from utils import combination_to_tuple


def s_basis_descriptors(homological_degree, n_vertices):
    """
    Return S-basis element descriptors at a given homological degree.

    These are the "labels" that identify each summand of the resolution:
      - degree 0: [()]
      - degree h >= 1: [((x_deg, y_deg), vertex_tuple), ...]

    Args:
        homological_degree: Homological degree h
        n_vertices: Number of vertices n in K_n

    Returns:
        List of S-basis descriptors
    """
    if homological_degree == 0:
        return [()]
    return [
        ((xd, homological_degree + 1 - xd), tuple(vc))
        for xd in range(1, homological_degree + 1)
        for vc in Combinations(range(n_vertices), homological_degree + 1)
    ]


def compute_basis_elements(homological_degree, ring_degree, n_vertices):
    """
    Compute k-basis elements for given homological and ring degrees.

    Each basis element is a pair (monomial_coeff, s_basis_descriptor).

    Args:
        homological_degree: Homological degree in the resolution
        ring_degree: Ring degree of the elements
        n_vertices: Number of vertices in K_n

    Returns:
        List of basis elements as tuples (monomial_coeff, s_basis_descriptor)
    """
    coeff_degree = ring_degree if homological_degree == 0 else ring_degree - (homological_degree + 1)
    if coeff_degree < 0:
        return []

    s_basis = s_basis_descriptors(homological_degree, n_vertices)

    coeff_combinations = Combinations(
        range(coeff_degree + 2 * n_vertices - 1), 2 * n_vertices - 1
    )

    return [
        (combination_to_tuple(comb, coeff_degree), s_elem)
        for comb in coeff_combinations
        for s_elem in s_basis
    ]
