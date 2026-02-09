"""
Basis construction for minimal free resolutions.

This module provides functions to construct basis elements for the minimal
free resolution of the binomial edge ideal of a complete graph.
"""

from sage.all import Combinations


def combination_to_tuple(combination, total_degree):
    """
    Convert a combination to a tuple representing monomial degrees.
    
    Args:
        combination: A tuple of integers representing positions
        total_degree: The total degree sum required
    
    Returns:
        A tuple with length len(combination)+1 and entries summing to total_degree
    """
    if not combination:
        return (total_degree,)
    
    t = (combination[0],)
    j = combination[0]
    for i in combination[1:]:
        t = t + (i - j - 1,)
        j = i
    t = t + (total_degree + len(combination) - (i + 1),)
    
    assert len(t) == len(combination) + 1
    assert sum(t) == total_degree
    
    return t


def compute_basis_elements(homological_degree, ring_degree, n_vertices):
    """
    Compute k-basis elements for given homological and ring degrees.
    
    Args:
        homological_degree: Homological degree in the resolution
        ring_degree: Ring degree of the elements
        n_vertices: Number of vertices in the complete graph
    
    Returns:
        List of basis elements as tuples (monomial_coeff, s_basis_element)
    """
    # Coefficient degree from the polynomial ring S
    coeff_degree = ring_degree if homological_degree == 0 else ring_degree - (homological_degree + 1)
    
    if coeff_degree < 0:
        return []
    
    # S-basis elements
    if homological_degree == 0:
        s_basis = [()]
    else:
        vertex_combinations = Combinations(range(n_vertices), homological_degree + 1)
        s_basis = [
            ((x_deg, homological_degree + 1 - x_deg), tuple(vertex_comb))
            for x_deg in range(1, homological_degree + 1)
            for vertex_comb in vertex_combinations
        ]
    
    # Coefficient monomials
    coeff_combinations = Combinations(range(coeff_degree + 2*n_vertices - 1), 2*n_vertices - 1)
    
    basis = []
    for comb in coeff_combinations:
        for s_elem in s_basis:
            # Create tuple of length 2*n_vertices
            coeff_tuple = combination_to_tuple(comb, coeff_degree)
            basis.append((coeff_tuple, s_elem))
    
    return basis


def compute_s_basis_elements(homological_degree, n_vertices):
    """
    Compute S-basis elements for a given homological degree.
    
    Args:
        homological_degree: Homological degree
        n_vertices: Number of vertices
    
    Returns:
        List of S-basis elements
    """
    if homological_degree == 0:
        return [(2*n_vertices*(0,), ())]
    
    vertex_combinations = Combinations(range(n_vertices), homological_degree + 1)
    return [
        (2*n_vertices*(0,), ((x_deg, homological_degree + 1 - x_deg), tuple(vertex_comb)))
        for x_deg in range(1, homological_degree + 1)
        for vertex_comb in vertex_combinations
    ]
