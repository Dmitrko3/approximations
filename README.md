# Approximations
5 Functions of Approximations


## Overview

This project was developed as part of the Numerical Analysis course in the second year of the Software Engineering degree.

The goal of the project is to implement and compare several interpolation and approximation methods used for estimating unknown function values based on a given set of known points.

Each method receives:

* A table of known points ((X,Y))
* A target value (X) that does not necessarily appear in the table

The algorithms then calculate an approximation of the corresponding function value using the selected interpolation method.

---

# Implemented Methods

## 1. Lagrange Interpolation

Implements the classical Lagrange interpolation formula.
The method builds a polynomial that passes exactly through all the given points using Lagrange basis polynomials.

### Main Features

* Accurate polynomial interpolation
* Uses all given points
* Includes validation for duplicate X values
* Returns exact values for known table points

---

## 2. Neville Interpolation

Implements Neville’s recursive interpolation algorithm.
The method gradually builds interpolation values using recursive combinations of nearby points until reaching the final approximation.

### Main Features

* Recursive interpolation process
* Numerically stable for many interpolation cases
* Efficient step-by-step computation
* Suitable for evaluating a single interpolation point

---

## 3. Linear Interpolation

Implements linear interpolation between two nearby points.
The method estimates the unknown value using a straight line connecting the closest surrounding points.

### Main Features

* Simple and fast computation
* Uses the nearest surrounding points only
* Handles extrapolation cases when needed
* Suitable for small local approximations

---

## 4. Polynomial Interpolation

Implements polynomial interpolation for estimating function values using a polynomial constructed from all known points.

### Main Features

* Uses all available data points
* Produces smooth approximations
* Accurate for polynomial-based datasets
* Includes input validation and error handling

---

## 5. Spline Interpolation

Implements spline interpolation using piecewise polynomials.
The method creates smooth transitions between points while preserving continuity.

### Main Features

* Smooth interpolation curves
* Better stability for larger datasets
* Reduces oscillation problems
* Produces continuous transitions between intervals

---

# Project Structure

```text
project/
│
├── lagrange.py
├── neville.py
├── linear.py
├── polynomial.py
├── spline.py
└── README.md
```

Each file contains:

* A complete implementation of the interpolation method
* Input validation
* Clear documentation and comments
* Example executions and output demonstrations

---

# Academic Purpose

This project was created for educational purposes as part of the Numerical Analysis course.

The implementations were designed to:

* Demonstrate understanding of interpolation algorithms
* Practice mathematical programming in Python
* Strengthen problem-solving and software development skills
* Apply numerical analysis concepts in practical coding tasks


