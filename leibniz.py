"""
Leibniz rule computations for DGA structure.

This module provides functions to compute the graded Leibniz rule
which must be satisfied by a DGA structure.
"""


def element_degree(element):
    """
    Compute the total degree of an element.
    
    Args:
        element: Element in the free module
    
    Returns:
        Total degree (x_deg + y_deg)
    """
    # Get the first monomial in the element
    if element.is_zero():
        return 0
    
    first_term = tuple(element)[0][0]
    
    if first_term[1] == ():
        return 0
    
    (x_deg, y_deg), _ = first_term[1]
    return x_deg + y_deg


def compute_leibniz_expression(elem1, elem2, differential, product, tensor_product):
    """
    Compute the graded Leibniz expression for two elements.
    
    The Leibniz rule states:
    d(f * g) = d(f) * g ± f * d(g)
    
    where the sign is (-1)^|f| for graded commutativity.
    
    Args:
        elem1: First element
        elem2: Second element
        differential: Differential morphism
        product: Product morphism
        tensor_product: Tensor product function
    
    Returns:
        d(elem1) * elem2 ± elem1 * d(elem2)
    """
    deg_f = element_degree(elem1)
    
    # Compute d(f) ⊗ g
    df_tensor_g = tensor_product([differential(elem1), elem2])
    
    # Compute f ⊗ d(g)
    f_tensor_dg = tensor_product([elem1, differential(elem2)])
    
    # Apply product
    product_df_g = product(df_tensor_g)
    product_f_dg = product(f_tensor_dg)
    
    # Graded Leibniz rule with sign
    if deg_f % 2 == 0:
        leibniz_expr = product_df_g - product_f_dg
    else:
        leibniz_expr = product_df_g + product_f_dg
    
    return leibniz_expr


def verify_leibniz_rule(elem1, elem2, differential, product, tensor_product):
    """
    Verify if the Leibniz rule holds for two elements.
    
    Args:
        elem1: First element
        elem2: Second element  
        differential: Differential morphism
        product: Product morphism
        tensor_product: Tensor product function
    
    Returns:
        True if Leibniz rule is satisfied, False otherwise
    """
    leibniz_expr = compute_leibniz_expression(
        elem1, elem2, differential, product, tensor_product
    )
    
    # Also compute d(product(elem1, elem2))
    elem1_tensor_elem2 = tensor_product([elem1, elem2])
    product_result = product(elem1_tensor_elem2)
    d_product = differential(product_result)
    
    # Check if they're equal
    difference = d_product - leibniz_expr
    
    return difference.is_zero()
