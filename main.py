"""
Demonstration of DGA resolution computations for binomial edge ideals.

This script walks through the core functionality:
  1. Building the free module and its basis
  2. Computing the differential
  3. Computing products with unknown coefficients
  4. Applying the Leibniz rule to derive equations
  5. Solving for unknown coefficients
"""

from sage.all import PolynomialRing, CombinatorialFreeModule, QQ, tensor
from basis import compute_basis_elements, s_basis_descriptors
from differential import create_differential_morphism
from product import create_product_morphism
from leibniz import compute_leibniz_expression
from solver import extract_equations, solve_equations_char0

def main():
    # Configuration

    n = 5          # vertices in K_n
    ring_deg = 5   # ring degree for basis elements

    # Polynomial ring with A, B, C variable families
    var_names = (
        [f'A{i}' for i in range(50)]
        + [f'B{i}' for i in range(50)]
        + [f'C{i}' for i in range(10)]
    )
    R = PolynomialRing(QQ, var_names)
    gens = R.gens()
    A_vars = gens[:50]
    B_vars = gens[50:100]
    C_vars = gens[100:110]

# Build the free module

    print("=" * 60)
    print("Building the free module F")
    print("=" * 60)

    all_basis = []
    for hdeg in range(n):
        elems = compute_basis_elements(hdeg, ring_deg, n)
        all_basis.extend(elems)
        print(f"Homological degree {hdeg}: {len(elems)} basis elements")

    F = CombinatorialFreeModule(R, tuple(all_basis))
    print(f"Total basis elements: {len(all_basis)}")

# S-basis descriptors

    print("\n" + "=" * 60)
    print("S-basis descriptors")
    print("=" * 60)

    for hdeg in range(min(3, n)):
        descs = s_basis_descriptors(hdeg, n)
        print(f"Degree {hdeg}: {len(descs)} descriptors")
        for d in descs[:5]:
            print(f"{d}")
        if len(descs) > 5:
            print(f"... ({len(descs) - 5} more)")

# The differential

    print("\n" + "=" * 60)
    print("Computing the differential")
    print("=" * 60)

    d = create_differential_morphism(F, n)
    zero_mon = (0,) * (2 * n)

    # Degree 1 example: d(e_{1,2} f_{1,1}) = x_1 y_2 - x_2 y_1
    f_deg1 = F.monomial((zero_mon, ((1, 1), (0, 1))))
    print(f"f = f_{{1,1}} e_{{1,2}}")
    print(f"d(f) = {d(f_deg1)}")

    # Degree 2 example
    f_deg2 = F.monomial((zero_mon, ((1, 2), (0, 1, 2))))
    print(f"g = f_{{1,2}} e_{{1,2,3}}")
    print(f"d(g) = {d(f_deg2)}")

    # Verify d^2 = 0
    print(f"d^2(g) = {d(d(f_deg2))}")
    verify_leibniz = d(d(f_deg2)).is_zero()
    print(f"d^2(g) = 0 verified: {verify_leibniz}")

# Multigraded product with unknown coefficients

    print("\n" + "=" * 60)
    print("Multigraded product with unknown coefficients")
    print("=" * 60)

    FF = tensor([F, F])
    mul_unknown = create_product_morphism(FF, F, A_vars, B_vars, n)

    f1 = F.monomial((zero_mon, ((1, 1), (0, 1))))
    f2 = F.monomial((zero_mon, ((1, 1), (1, 2))))
    product_fg = mul_unknown(tensor([f1, f2]))

    print(f"f1 = f_{{1,1}} e_{{1,2}}")
    print(f"f2 = f_{{1,1}} e_{{2,3}}")
    print(f"f1 * f2 ={product_fg}")

# Leibniz rule and equation extraction

    print("\n" + "=" * 60)
    print("Leibniz rule: d(f*g) = d(f)*g + (-1)^|f| f*d(g)")
    print("=" * 60)

    d_product = d(product_fg)
    leibniz_expr = compute_leibniz_expression(
        f1, f2, d, mul_unknown, tensor
    )
    difference = d_product - leibniz_expr

    print(f"d(f1 * f2) - Leibniz(f1,f2) = {difference}")

    # Extract and solve
    equations = extract_equations(difference)
    print(f"Extracted {len(equations)} equations")

    solution = solve_equations_char0(equations)
    if solution is not None:
        print(f"Solution: {solution}")
    else:
        print("No solution found (or trivially satisfied)")

# Cross-degree product example

    print("\n" + "=" * 60)
    print("Cross-degree product (degree 1 * degree 2)")
    print("=" * 60)

    f1 = F.monomial((zero_mon, ((1, 1), (0, 1))))
    f2 = F.monomial((zero_mon, ((1, 2), (2, 3, 4))))

    product_cross = mul_unknown(tensor([f1, f2]))
    print(f"f1 = f_{{1,1}} e_{{1,2}} (degree 1)")
    print(f"f2 = f_{{1,2}} e_{{3,4,5}} (degree 2)")
    print(f"f1 * f2 (with unknowns) = {product_cross}")

    d_cross = d(product_cross)
    leibniz_cross = compute_leibniz_expression(
        f1, f2, d, mul_unknown, tensor
    )
    diff_cross = d_cross - leibniz_cross

    equations_cross = extract_equations(diff_cross)
    print(f"Extracted {len(equations_cross)} equations from Leibniz rule")


if __name__ == "__main__":
    main()