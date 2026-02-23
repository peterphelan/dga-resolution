"""
Shared utilities for the DGA resolution project.

Collects helper functions used across multiple modules.
"""


def increment_tuple(tup, position):
    """Increment a tuple at a specific position by 1."""
    return tuple(
        tup[i] + 1 if i == position else tup[i]
        for i in range(len(tup))
    )


def add_tuples(t1, t2):
    """Componentwise addition of two tuples of the same length."""
    assert len(t1) == len(t2)
    return tuple(a + b for a, b in zip(t1, t2))


def scalar_multiply_basis(monomial_tuple, basis_element, free_module):
    """
    Multiply a basis element by a monomial (scalar multiplication in S).

    Args:
        monomial_tuple: Tuple of degrees for 2*n_vertices variables
        basis_element: Basis element (monomial_coeff, s_basis_part)
        free_module: The combinatorial free module

    Returns:
        Product as an element of the free module
    """
    new_monomial = add_tuples(monomial_tuple, basis_element[0])
    return free_module.monomial((new_monomial, basis_element[1]))


def vertex_overlap(e1_info, e2_info):
    """
    Count shared vertices between two S-basis descriptors.

    Products with overlap > 1 are always zero by multigrading, so
    they can be skipped in the computation loop.
    """
    if e1_info == () or e2_info == ():
        return 0
    return len(set(e1_info[1]) & set(e2_info[1]))


def combination_to_tuple(combination, total_degree):
    """
    Convert a combination to a tuple representing monomial degrees.

    Given a combination (subset of {0,...,total_degree + len(combination) - 1}),
    produce the tuple of gap lengths, which sums to total_degree.
    """
    if not combination:
        return (total_degree,)

    parts = [combination[0]]
    for i in range(1, len(combination)):
        parts.append(combination[i] - combination[i - 1] - 1)
    parts.append(total_degree + len(combination) - (combination[-1] + 1))

    assert len(parts) == len(combination) + 1
    assert sum(parts) == total_degree
    return tuple(parts)
