"""
Coefficient Manager for minimal free resolution.

Stores and retrieves the coefficients of terms in products of basis elements
determined by the code. Coefficients are persisted to JSON and reconstructed 
as polynomial ring elements when loaded.
"""

import json
from pathlib import Path
from typing import Dict, Tuple, Optional


class CoefficientManager:
    """
    Manages product coefficients for proposed DGA structures.

    Keys: pairs of S-basis element descriptors
          ((xdeg1, ydeg1), vertices1), ((xdeg2, ydeg2), vertices2)
    Values: tuples of polynomial ring elements (or strings if no ring).
    """

    def __init__(self, filepath=None, polynomial_ring=None):
        self.filepath = Path(filepath) if filepath else None
        self.polynomial_ring = polynomial_ring
        self.coefficients: Dict[Tuple, Tuple] = {}

        if self.filepath and self.filepath.exists():
            self.load()

    # Storage operations

    def add_coefficient(self, elem1_info, elem2_info, coeffs):
        """Store product coefficients."""
        key = (elem1_info, elem2_info)
        self.coefficients[key] = self._normalize_coeffs(coeffs)

    def get_coefficient(self, elem1_info, elem2_info):
        """Retrieve coefficients, checking both key orderings."""
        for key in [(elem1_info, elem2_info), (elem2_info, elem1_info)]:
            if key in self.coefficients:
                return self.coefficients[key]
        return None

    def has_coefficient(self, elem1_info, elem2_info):
        return self.get_coefficient(elem1_info, elem2_info) is not None

    def get_all_coefficients(self):
        return self.coefficients.copy()

    def clear(self):
        self.coefficients.clear()

    # Persistence

    def save(self, filepath=None):
        """Save to JSON. Ring elements are serialized as strings."""
        path = Path(filepath) if filepath else self.filepath
        if not path:
            raise ValueError("No filepath specified for saving")

        data = {
            json.dumps(k): [str(c) for c in v]
            for k, v in self.coefficients.items()
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)

    def load(self, filepath=None):
        """Load from JSON. Converts strings back to ring elements if ring provided."""
        path = Path(filepath) if filepath else self.filepath
        if not path or not path.exists():
            return

        with open(path, 'r') as f:
            data = json.load(f)

        to_ring = self.polynomial_ring or str
        self.coefficients = {
            self._lists_to_tuples(json.loads(k)): tuple(to_ring(c) for c in v)
            for k, v in data.items()
        }

    # Internal helpers

    def _normalize_coeffs(self, coeffs):
        """Unpack SymPy FiniteSet/Tuple and convert to ring elements."""
        from sympy import FiniteSet
        from sympy.core.containers import Tuple as SympyTuple

        if isinstance(coeffs, FiniteSet):
            coeffs = next(iter(coeffs))
        if isinstance(coeffs, SympyTuple):
            coeffs = tuple(coeffs)

        convert = self.polynomial_ring or str
        return tuple(convert(str(c)) for c in coeffs)

    @staticmethod
    def _lists_to_tuples(obj):
        """Recursively convert nested lists to tuples (for JSON key reconstruction)."""
        if isinstance(obj, list):
            return tuple(CoefficientManager._lists_to_tuples(x) for x in obj)
        return obj

    def __len__(self):
        return len(self.coefficients)

    def __repr__(self):
        return f"CoefficientManager(entries={len(self)}, filepath={self.filepath})"
