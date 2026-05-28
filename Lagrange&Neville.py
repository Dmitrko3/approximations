"""
Interpolation assignment:
- Lagrange interpolation
- Neville interpolation
"""

from typing import List, Tuple

Point = Tuple[float, float]


def validate_input(points: List[Point], target_x: float) -> float | None:
    """
    Validates that the interpolation input is legal.
    If target_x is already in the table, returns the known y-value.
    Otherwise returns None if validation passes.
    """
    # Check if target_x is already in the table
    for x, y in points:
        if x == target_x:
            return y

    if len(points) < 2:
        raise ValueError("At least two points are required.")

    x_values = [point[0] for point in points]

    if len(x_values) != len(set(x_values)):
        raise ValueError("X values must be unique.")

    return None
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
    validate_input(points, target_x)
    n = len(points)
    p = [point[1] for point in points] # מערך חד ממדי שמאתחל עם ערכי ה-Y

    for j in range(1, n):
        for i in range(n - j):
            xi = points[i][0]
            xj = points[i + j][0]
            p[i] = ((target_x - xj) * p[i] - (target_x - xi) * p[i + 1]) / (xi - xj)

    return p[0]


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


if __name__ == "__main__":
    main()
