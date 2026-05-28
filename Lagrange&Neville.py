"""
Interpolation assignment:
- Lagrange interpolation
- Neville interpolation
"""

from typing import List, Tuple

Point = Tuple[float, float]


def validate_input(points: List[Point], target_x: float) -> None:
    """
    Validates that the interpolation input is legal.
    """
    if len(points) < 2:
        raise ValueError("At least two points are required.")

    x_values = [point[0] for point in points]

    if len(x_values) != len(set(x_values)):
        raise ValueError("X values must be unique.")

    if target_x in x_values:
        raise ValueError("Target x should not already appear in the table.")


def lagrange_interpolation(points: List[Point], target_x: float) -> float:
    """
    Calculates interpolation value using the Lagrange method.
    """
    validate_input(points, target_x)

    result = 0.0

    for i in range(len(points)):
        xi, yi = points[i]
        basis = 1.0

        for j in range(len(points)):
            if i != j:
                xj, _ = points[j]
                basis *= (target_x - xj) / (xi - xj)

        result += yi * basis

    return result


def neville_interpolation(points: List[Point], target_x: float) -> float:
    """
    Calculates interpolation value using Neville's method.
    """
    validate_input(points, target_x)

    n = len(points)
    table = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        table[i][0] = points[i][1]

    for j in range(1, n):
        for i in range(n - j):
            xi = points[i][0]
            xj = points[i + j][0]

            table[i][j] = (
                (target_x - xj) * table[i][j - 1]
                - (target_x - xi) * table[i + 1][j - 1]
            ) / (xi - xj)

    return table[0][n - 1]


def main() -> None:
    """
    Main program:
    defines table points, defines target x,
    and prints the result of both interpolation methods.
    """

    points = [
        (1.0, 1.0),
        (2.0, 4.0),
        (3.0, 9.0),
        (4.0, 16.0),
    ]

    target_x = 2.5

    lagrange_result = lagrange_interpolation(points, target_x)
    neville_result = neville_interpolation(points, target_x)

    print("Interpolation Results")
    print("---------------------")
    print(f"Table points: {points}")
    print(f"Target x: {target_x}")
    print(f"Lagrange result: {lagrange_result}")
    print(f"Neville result: {neville_result}")


if __name__ == "main":
    main()
