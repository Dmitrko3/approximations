"""
Numerical Analysis - Assignment Runner
Compares Neville and Cubic Spline interpolation over 200 uniform points.
"""

from lagrange_neville import neville_interpolation
from cubic_spline import cubic_spline


def main():
    # The 8 original table points from the domain [0, 3]
    points = [
        (0.0, 0.500),
        (0.3, 0.779),
        (0.7, -0.054),
        (1.0, -0.591),
        (1.5, -0.368),
        (2.0, -0.715),
        (2.5, -0.056),
        (3.0, 0.976)
    ]

    # Split into x and y arrays for the cubic_spline function signature
    x_nodes = [p[0] for p in points]
    y_nodes = [p[1] for p in points]

    # Generate 200 uniformly distributed points in the domain [0, 3]
    num_points = 200
    start_x, end_x = 0.0, 3.0

    # Using (num_points - 1) ensures the last point is exactly 3.0
    test_points = [start_x + i * (end_x - start_x) / (num_points - 1) for i in range(num_points)]

    differences = []

    # 1. Calculate the difference at each test point
    for x in test_points:
        y_neville = neville_interpolation(points, x)
        y_spline = cubic_spline(x_nodes, y_nodes, x)

        # Absolute gap between the two methods
        diff = abs(y_neville - y_spline)
        differences.append(diff)

    # 2. Find the maximum difference (delta_max)
    delta_max = max(differences)

    # 3. Calculate the range of values (R) from the original 8 table values
    # R = max(f) - min(f)
    r_val = max(y_nodes) - min(y_nodes)

    # 4. Calculate delta_max / R
    ratio = delta_max / r_val

    # Print the assignment output
    print("=" * 40)
    print(" Interpolation Comparison Results")
    print("=" * 40)
    print(f"Test Points Evaluated : {num_points}")
    print(f"Max Difference (delta_max) : {delta_max:.6f}")
    print(f"Function Range (R)         : {r_val:.6f}")
    print(f"Ratio (delta_max / R)      : {ratio:.6f}")
    print("=" * 40)


if __name__ == "__main__":
    main()