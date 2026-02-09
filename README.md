# DGA Resolutions

This repository contains code for investigating the existence of multigraded DGA (Differential Graded Algebra) structures on the minimal free resolutions of binomial edge ideals of complete graphs.

## Current Features

- Construction of basis elements on minimal free resolution.
- Differential defined on the minimal free resolution.
- Multigraded product of basis elements with unknown coefficients.
- SageMath + Docker reproducible environment.

## Planned Features

- Extract systems of equations in unknown coefficients arising from the Leibniz rule.
- Reduce systems of equations in unknown coefficients.
- Coefficient manager to save, store, and load reduced systems of equations in JSON format.
- Further reduce systems of equations using associativity.
- Solve systems of equations, or verify non-existence of solutions.

## Usage

To deploy the Docker container, run:

```bash
docker build -t dga-resolution .
docker run dga-resolution
```
