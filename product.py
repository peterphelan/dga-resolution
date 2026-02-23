"""
Product structure for the minimal free resolution with unknown coefficients.

This module implements the multigraded product structure on the resolution,
with unknown coefficients represented as polynomial ring variables.
"""

from utils import scalar_multiply_basis, add_tuples


def _coeff_increment_helpers(elem1, elem2, n_vertices):
    """
    Build helpers that produce the combined monomial coefficient
    with one variable incremented (x_i or y_i contribution).
    """
    base = add_tuples(elem1[0], elem2[0])

    def with_x(i):
        return tuple(base[j] + (1 if j == i else 0) for j in range(len(base)))

    def with_y(i):
        return tuple(base[j] + (1 if j == n_vertices + i else 0) for j in range(len(base)))

    return with_x, with_y


def compute_multigraded_product(basis_pair, free_module, A_vars, B_vars, n_vertices):
    """
    Compute the multigraded product of two basis elements with unknown coefficients.

    Returns all terms permitted by the multigrading, with A/B unknowns.

    Args:
        basis_pair: Tuple (basis_elem1, basis_elem2)
        free_module: The free module
        A_vars: List of 'A' coefficient variables
        B_vars: List of 'B' coefficient variables
        n_vertices: Number of vertices
    """
    elem1, elem2 = basis_pair

    # Scalar multiplication when one factor is in degree 0
    if elem1[1] == ():
        return scalar_multiply_basis(elem1[0], elem2, free_module)
    if elem2[1] == ():
        return scalar_multiply_basis(elem2[0], elem1, free_module)

    # Check if product exceeds maximum homological degree
    hdeg1 = len(elem1[1][1])
    hdeg2 = len(elem2[1][1])
    if hdeg1 + hdeg2 > n_vertices + 1:
        return free_module.zero()

    (x1, y1), vdegs1 = elem1[1]
    (x2, y2), vdegs2 = elem2[1]

    # Graded commutativity: ensure canonical ordering
    if hdeg1 == 2 and hdeg2 == 2 and vdegs1[0] > vdegs2[0]:
        return -compute_multigraded_product(
            (elem2, elem1), free_module, A_vars, B_vars, n_vertices
        )

    new_vdegs = tuple(sorted(set(vdegs1) | set(vdegs2)))
    new_x = x1 + x2
    new_y = y1 + y2
    intersection = set(vdegs1) & set(vdegs2)
    coeff_x, coeff_y = _coeff_increment_helpers(elem1, elem2, n_vertices)

    result = free_module.zero()

    if len(intersection) == 0:
        for k, vertex in enumerate(new_vdegs):
            dropped = new_vdegs[:k] + new_vdegs[k + 1:]
            result += A_vars[k + 1] * free_module.monomial(
                (coeff_x(vertex), ((new_x - 1, new_y), dropped))
            )
            result += B_vars[k + 1] * free_module.monomial(
                (coeff_y(vertex), ((new_x, new_y - 1), dropped))
            )

    elif len(intersection) == 1:
        for vertex in intersection:
            k = new_vdegs.index(vertex)
            result += A_vars[k + 1] * free_module.monomial(
                (coeff_x(vertex), ((new_x - 1, new_y), new_vdegs))
            )
            result += B_vars[k + 1] * free_module.monomial(
                (coeff_y(vertex), ((new_x, new_y - 1), new_vdegs))
            )

    # len(intersection) > 1: zero by multigrading

    return result


def create_product_morphism(tensor_module, free_module, A_vars, B_vars, n_vertices):
    """Create the symbolic product as a module morphism F ⊗ F → F."""
    def product_on_basis(basis_pair):
        return compute_multigraded_product(
            basis_pair, free_module, A_vars, B_vars, n_vertices
        )
    return tensor_module.module_morphism(product_on_basis, codomain=free_module)
