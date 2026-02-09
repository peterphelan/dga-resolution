"""
Differential map for minimal free resolutions.

This module implements the differential on the minimal free resolution
of the binomial edge ideal of a complete graph.
"""


def increment_tuple(tup, position):
    """
    Increment a tuple at a specific position.
    
    Args:
        tup: Input tuple
        position: Index to increment
    
    Returns:
        New tuple with value at position incremented by 1
    """
    return tuple(
        tup[i] + 1 if i == position else tup[i]
        for i in range(len(tup))
    )


def compute_differential(basis_element, free_module, n_vertices):
    """
    Compute the differential of a basis element.
    
    Args:
        basis_element: Tuple (monomial_coeff, s_basis_part)
        free_module: The combinatorial free module
        n_vertices: Number of vertices
    
    Returns:
        Element in the free module representing d(basis_element)
    """
    monomial_coeff, s_basis_part = basis_element
    
    # Homological degree 0: differential is zero
    if s_basis_part == ():
        return free_module.zero()
    
    (x_deg, y_deg), vertex_degs = s_basis_part
    
    # Homological degree 1: special case
    if len(vertex_degs) == 2:
        i, j = vertex_degs
        term1 = free_module.monomial(
            (increment_tuple(increment_tuple(monomial_coeff, i), j + n_vertices), ())
        )
        term2 = free_module.monomial(
            (increment_tuple(increment_tuple(monomial_coeff, j), i + n_vertices), ())
        )
        return term1 - term2
    
    # Homological degree > 1: general formula
    result = free_module.zero()
    
    for idx, vertex in enumerate(vertex_degs):
        sign_x = 1 if idx % 2 == 0 else -1
        sign_y = -1 if idx % 2 == 0 else 1
        
        remaining_vertices = vertex_degs[:idx] + vertex_degs[idx+1:]
        
        if x_deg > 1:
            new_monomial = increment_tuple(monomial_coeff, vertex)
            new_s_basis = ((x_deg - 1, y_deg), remaining_vertices)
            result += sign_x * free_module.monomial((new_monomial, new_s_basis))
        
        if y_deg > 1:
            new_monomial = increment_tuple(monomial_coeff, vertex + n_vertices)
            new_s_basis = ((x_deg, y_deg - 1), remaining_vertices)
            result += sign_y * free_module.monomial((new_monomial, new_s_basis))
    
    return result


def create_differential_morphism(free_module, n_vertices):
    """
    Create the differential as a module morphism.
    
    Args:
        free_module: The combinatorial free module
        n_vertices: Number of vertices
    
    Returns:
        Module morphism representing the differential
    """
    def diff_on_basis(b):
        return compute_differential(b, free_module, n_vertices)
    
    return free_module.module_morphism(diff_on_basis, codomain=free_module)
