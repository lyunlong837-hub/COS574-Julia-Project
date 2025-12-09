# Implementation and Evaluation of Nonlinear Equation Solvers by Julia

This project implements several classic **nonlinear equation solvers** in Julia
and provides a simple command-line **calculator** to find roots of equations

\[
f(x) = 0
\]

Methods:
- **Bisection method**
- **Newton method**
- **Secant method**
- **Broyden method (1D quasi-Newton)**

---

## 1. Requirements

- **Julia** ≥ 1.8 (recommended 1.10+)
- No external packages are required.  
  The project only uses Julia Base.

---

## 2. Project Structure

```text
.
├── nonlinear_calculator.jl   # Main code: solvers + CLI "calculator"
└── README.md                 # This file
```

## 3. How to Run the Project

In the project directory:
```
julia nonlinear_calculator.jl
```

Input the function f(x) as a Julia expression, for example:
```
x^3 - x - 2

sin(x) - 0.5

exp(x) - 3x
```

Choose method:
```
1 → Bisection

2 → Newton

3 → Secant

4 → Broyden (1D)
```

Provide parameters:
```
Tolerance (e.g. 1e-8)

Maximum iterations (e.g. 100)

Initial interval or initial guesses, depending on the method
```

Get result:
```
Approximate root

f(root)

Number of iterations

Whether the method converged
```