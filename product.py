"""
Product structure for the DGA with unknown coefficients.

This module implements the multigraded product structure on the resolution,
with unknown coefficients represented as polynomial ring variables.
"""


def scalar_multiply_basis(monomial_tuple, basis_element, free_module):
    """
    Multiply a basis element by a monomial (scalar multiplication).
    
    Args:
        monomial_tuple: Tuple of degrees for 2*n_vertices variables
        basis_element: Basis element to multiply
        free_module: The free module
    
    Returns:
        Product as an element of the free module
    """
    assert len(monomial_tuple) == len(basis_element[0])
    
    new_monomial = tuple(
        monomial_tuple[i] + basis_element[0][i]
        for i in range(len(monomial_tuple))
    )
    
    return free_module.monomial((new_monomial, basis_element[1]))


def compute_multigraded_product(basis_pair, free_module, coeff_vars_a, coeff_vars_b, n_vertices):
    """
    Compute the multigraded product of two basis elements with unknown coefficients.
    
    This function returns all terms permitted by the multigrading with unknown
    coefficients from the polynomial ring.
    
    Args:
        basis_pair: Tuple (basis_elem1, basis_elem2)
        free_module: The free module
        coeff_vars_a: List of 'A' coefficient variables
        coeff_vars_b: List of 'B' coefficient variables
        n_vertices: Number of vertices
    
    Returns:
        Element in the free module with coefficient variables
    """
    elem1, elem2 = basis_pair
    
    # Base case: one factor has homological degree 0
    if elem1[1] == ():
        return scalar_multiply_basis(elem1[0], elem2, free_module)
    if elem2[1] == ():
        return scalar_multiply_basis(elem2[0], elem1, free_module)
    
    # Check if product exceeds maximum homological degree
    hdeg1 = len(elem1[1][1])
    hdeg2 = len(elem2[1][1])
    if hdeg1 + hdeg2 > n_vertices + 1:
        return free_module.zero()
    
    result = free_module.zero()
    
    # Extract degree information
    (x_deg1, y_deg1), vdegs1 = elem1[1]
    (x_deg2, y_deg2), vdegs2 = elem2[1]
    
    # Graded commutativity: ensure ordering
    if hdeg1 == 2 and hdeg2 == 2 and vdegs1[0] > vdegs2[0]:
        return -compute_multigraded_product(
            (elem2, elem1), free_module, coeff_vars_a, coeff_vars_b, n_vertices
        )
    
    # Combined vertex degrees and intersections
    new_vdegs = tuple(sorted(set(vdegs1) | set(vdegs2)))
    new_x_deg = x_deg1 + x_deg2
    new_y_deg = y_deg1 + y_deg2
    intersection = set(vdegs1) & set(vdegs2)
    
    # Helper functions for new coefficients
    def coeff_with_x_increment(i):
        return tuple(
            elem1[0][j] + elem2[0][j] + (1 if j == i else 0)
            for j in range(len(elem1[0]))
        )
    
    def coeff_with_y_increment(i):
        return tuple(
            elem1[0][j] + elem2[0][j] + (1 if j == n_vertices + i else 0)
            for j in range(len(elem1[0]))
        )
    
    # Case: no vertex intersection
    if len(intersection) == 0:
        for k, vertex in enumerate(new_vdegs):
            modified_vdegs = new_vdegs[:k] + new_vdegs[k+1:]
            
            coeff_x = coeff_with_x_increment(vertex)
            coeff_y = coeff_with_y_increment(vertex)
            
            result += coeff_vars_a[k+1] * free_module.monomial(
                (coeff_x, ((new_x_deg - 1, new_y_deg), modified_vdegs))
            )
            result += coeff_vars_b[k+1] * free_module.monomial(
                (coeff_y, ((new_x_deg, new_y_deg - 1), modified_vdegs))
            )
    
    # Case: exactly one vertex intersection
    elif len(intersection) == 1:
        for vertex in intersection:
            k = new_vdegs.index(vertex)
            
            coeff_x = coeff_with_x_increment(vertex)
            coeff_y = coeff_with_y_increment(vertex)
            
            result += coeff_vars_a[k+1] * free_module.monomial(
                (coeff_x, ((new_x_deg - 1, new_y_deg), new_vdegs))
            )
            result += coeff_vars_b[k+1] * free_module.monomial(
                (coeff_y, ((new_x_deg, new_y_deg - 1), new_vdegs))
            )
    
    # Case: more than one vertex intersection (violates multigrading)
    else:
        result = free_module.zero()
    
    return result


def create_product_morphism(tensor_module, free_module, coeff_vars_a, coeff_vars_b, n_vertices):
    """
    Create the product as a module morphism from F ⊗ F → F.
    
    Args:
        tensor_module: Tensor product module F ⊗ F
        free_module: Target free module F
        coeff_vars_a: List of 'A' coefficient variables
        coeff_vars_b: List of 'B' coefficient variables
        n_vertices: Number of vertices
    
    Returns:
        Module morphism representing the product
    """
    def product_on_basis(basis_pair):
        return compute_multigraded_product(
            basis_pair, free_module, coeff_vars_a, coeff_vars_b, n_vertices
        )
    
    return tensor_module.module_morphism(product_on_basis, codomain=free_module)
