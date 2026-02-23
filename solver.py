"""
Equation extraction and solving for Leibniz constraints.

Extracts the linear system from a free-module expression and solves it,
substituting C-variables for any remaining free parameters.
"""

import re
from sympy import simplify, Eq, linsolve, linear_eq_to_matrix, Symbol, FiniteSet
from sympy.core.containers import Tuple as SympyTuple


def extract_equations(expr):
    """
    Extract coefficient equations from a free module expression.

    Each monomial in the free module contributes one equation (its
    coefficient must vanish).

    Args:
        expr: Free module element whose coefficients should all be zero

    Returns:
        List of equation strings
    """
    equations = []
    expr_str = str(expr)
    pattern = r'\(([^()]+)\)'
    for m in re.findall(pattern, expr_str):
        if re.match(r'^[A-Za-z0-9+\-*/^ ]+$', m):
            equations.append(m)
    return equations


def solve_equations_char0(equations):
    """
    Solve a linear system over Q using SymPy.

    Args:
        equations: List of expression strings that should equal zero

    Returns:
        SymPy FiniteSet solution, or None if no solution / empty system
    """
    if not equations:
        return None

    sympy_eqs = []
    for expr_str in equations:
        try:
            sympy_eqs.append(Eq(simplify(expr_str), 0))
        except Exception:
            continue

    if not sympy_eqs:
        return None

    syms = sorted(
        set().union(*(eq.free_symbols for eq in sympy_eqs)),
        key=str,
    )
    A, b = linear_eq_to_matrix(sympy_eqs, syms)
    sol = linsolve((A, b), syms)
    return sol if sol else None


def solve_equations_finite_char(equations, poly_ring):
    """
    Solve a linear system over a finite field using Sage's Gröbner bases.

    Args:
        equations: List of expression strings
        poly_ring: Sage polynomial ring over the finite field

    Returns:
        Gröbner basis of the ideal
    """
    polys = [poly_ring(eq) for eq in equations]
    return poly_ring.ideal(polys).groebner_basis()


def substitute_free_variables(solution, next_c_index):
    """
    Replace free B-variables in a solution with C-variables.

    C-variables are used to represent free parameters in a
    reduced system of equations.

    Args:
        solution: SymPy solution (FiniteSet or tuple)
        next_c_index: Starting index for new C-variables

    Returns:
        (new_solution, updated_next_c_index)
    """
    if solution is None:
        return None, next_c_index

    # Unwrap FiniteSet
    if isinstance(solution, FiniteSet):
        sol_tuple = next(iter(solution))
    elif isinstance(solution, SympyTuple):
        sol_tuple = solution
    else:
        return solution, next_c_index

    # Find free B-variables
    all_syms = set()
    for e in sol_tuple:
        all_syms |= e.free_symbols
    b_vars = sorted(
        [s for s in all_syms if str(s).startswith('B')],
        key=lambda x: int(str(x)[1:]),
    )

    if not b_vars:
        return solution, next_c_index

    subs = {}
    for bv in b_vars:
        subs[bv] = Symbol(f'C{next_c_index}')
        next_c_index += 1

    result = tuple(e.subs(subs) for e in sol_tuple)
    return FiniteSet(SympyTuple(*result)), next_c_index
