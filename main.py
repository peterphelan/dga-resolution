"""
Main script demonstrating DGA resolution computations, including:
- Computing multigraded products with unknown coefficients
- Applying the Leibniz rule
"""

from sage.all import PolynomialRing, CombinatorialFreeModule, QQ, GF, tensor
from product import create_product_morphism
from differential import create_differential_morphism
from leibniz import compute_leibniz_expression
from basis import compute_basis_elements


def main():
    """Main function demonstrating all functionality."""

    # Configuration
    n_vertices = 5
    r_degree = 5
    base_field = QQ
    
    # Create polynomial ring with A, B, C variables
    var_names = (
        [f'A{i}' for i in range(50)] +
        [f'B{i}' for i in range(50)] +
        [f'C{i}' for i in range(10)]
    )
    R = PolynomialRing(base_field, var_names)
    variables = R.gens()
    
    A_vars = variables[:50]
    B_vars = variables[50:100]
    C_vars = variables[100:110]
    
    # Build free module F
    all_basis = []
    for hdeg in range(n_vertices):
        basis_at_degree = compute_basis_elements(hdeg, r_degree, n_vertices)
        all_basis.extend(basis_at_degree)
    
    F = CombinatorialFreeModule(R, tuple(all_basis))
    FF = tensor([F, F])
    
    # Create morhism for differential
    differential = create_differential_morphism(F, n_vertices)

    # Create morphism for multigraded product with unknown coefficients
    product_unknown = create_product_morphism(FF, F, A_vars, B_vars, n_vertices)

    # Choose basis elements for demonstration
    zero_mon = (0,) * (2 * n_vertices)
    f1 = F.monomial((zero_mon, ((1, 1), (0, 1))))
    f2 = F.monomial((zero_mon, ((1, 2), (1, 2, 3))))
    print(f"First basis element: {f1}")
    print(f"Second basis element: {f2}")

    # Test differential by applying it to the basis elements
    print("Applying differential to basis elements...")
    df1 = differential(f1)
    df2 = differential(f2)
    print(f"d(f1) = {df1}")
    print(f"d(f2) = {df2}")
    
    # Compute product with unknowns
    print(f"Computing multigraded product of f1 and f2 with unknown coefficients...")
    product_fg = product_unknown(tensor([f1, f2]))
    print(f"Product is {product_fg}")
    
    # Compute the differential of the product and apply the Leibniz rule
    print("Applying Leibniz rule to the product...")
    d_product = differential(product_fg)
    leibniz_expr = compute_leibniz_expression(
        f1, f2, differential, product_unknown, tensor
    )
    print(f"d(f1 * f2) = {d_product}")
    

if __name__ == "__main__":
   main()
