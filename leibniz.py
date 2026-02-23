"""
Leibniz rule computations for minimal free resolution.

Provides functions to compute and verify the graded Leibniz rule:
    d(f * g) = d(f) * g + (-1)^|f| * f * d(g)
"""


def element_degree(element):
    """
    Compute the homological degree of an element.

    Returns x_deg + y_deg from the first term, or 0 for degree-0 elements.
    """
    if element.is_zero():
        return 0

    first_term = tuple(element)[0][0]

    if first_term[1] == ():
        return 0

    (x_deg, y_deg), _ = first_term[1]
    return x_deg + y_deg


def compute_leibniz_expression(elem1, elem2, differential, product, tensor_product):
    """
    Compute d(f)*g + (-1)^|f| f*d(g).

    This is the right-hand side of the Leibniz rule; by the definition 
    of a DGA, it must be equal to d(f*g).

    Args:
        elem1, elem2: Elements of the free module
        differential: Differential morphism d
        product: Product morphism (symbolic or concrete)
        tensor_product: Tensor product function
    """
    deg_f = element_degree(elem1)

    df_g = product(tensor_product([differential(elem1), elem2]))
    f_dg = product(tensor_product([elem1, differential(elem2)]))

    if deg_f % 2 == 0:
        return df_g - f_dg
    else:
        return df_g + f_dg


def verify_leibniz_rule(elem1, elem2, differential, product, tensor_product):
    """
    Verify d(f*g) = d(f)*g + (-1)^|f| * f*d(g) for two elements.

    Returns True if the Leibniz rule is satisfied.
    """
    leibniz_rhs = compute_leibniz_expression(
        elem1, elem2, differential, product, tensor_product
    )

    fg = product(tensor_product([elem1, elem2]))
    d_fg = differential(fg)

    return (d_fg - leibniz_rhs).is_zero()
