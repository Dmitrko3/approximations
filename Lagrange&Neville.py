"""
Interpolation assignment:
- Lagrange interpolation
- Neville interpolation
"""

from typing import List, Tuple

import matplotlib.pyplot as plt
from cubic_spline import cubic_spline

Point = Tuple[float, float]


def validate_input(points: List[Point], target_x: float) -> float | None:
    """
    Validates that the interpolation input is legal.
    If target_x is already in the table, returns the known y-value.
    Otherwise, returns None if validation passes.
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
    # """
    # Main program:
    # defines table points, defines target x,
    # and prints the result of both interpolation methods.
    # """
    #
    # points = [
    #     (0.0, 0.500),
    #     (0.3, 0.779),
    #     (0.7, -0.054),
    #     (1.0, -0.591),
    #     (1.5, -0.368),
    #     (2.0, -0.715),
    #     (2.5, -0.056),
    #     (3.0, 0.976)
    # ]
    #
    # target_x = 2.5
    #
    # lagrange_result = lagrange_interpolation(points, target_x)
    # neville_result = neville_interpolation(points, target_x)
    #
    # print("Interpolation Results")
    # print("---------------------")
    # print(f"Table points: {points}")
    # print(f"Target x: {target_x}")
    # #print(f"Lagrange result: {lagrange_result}")
    # print(f"Neville result: {neville_result}")

#     def main():
#         # The 8 original table points from the domain [0, 3]
#         points = [
#             (0.0, 0.500),
#             (0.3, 0.779),
#             (0.7, -0.054),
#             (1.0, -0.591),
#             (1.5, -0.368),
#             (2.0, -0.715),
#             (2.5, -0.056),
#             (3.0, 0.976)
#         ]
#
#         # Split into x and y arrays for the cubic_spline function signature
#         x_nodes = [p[0] for p in points]
#         y_nodes = [p[1] for p in points]
#
#         # Generate 200 uniformly distributed points in the domain [0, 3]
#         num_points = 200
#         start_x, end_x = 0.0, 3.0
#
#         # Using (num_points - 1) ensures the last point is exactly 3.0
#         test_points = [start_x + i * (end_x - start_x) / (num_points - 1) for i in range(num_points)]
#
#         differences = []
#
#         # 1. Calculate the difference at each test point
#         for x in test_points:
#             y_neville = neville_interpolation(points, x)
#             y_spline = cubic_spline(x_nodes, y_nodes, x)
#
#             # Absolute gap between the two methods
#             diff = abs(y_neville - y_spline)
#             differences.append(diff)
#
#         # 2. Find the maximum difference (delta_max)
#         delta_max = max(differences)
#
#         # 3. Calculate the range of values (R) from the original 8 table values
#         # R = max(f) - min(f)
#         r_val = max(y_nodes) - min(y_nodes)
#
#         # 4. Calculate delta_max / R
#         ratio = delta_max / r_val
#
#         # Print the assignment output
#         print("=" * 40)
#         print(" Interpolation Comparison Results")
#         print("=" * 40)
#         print(f"Test Points Evaluated : {num_points}")
#         print(f"Max Difference (delta_max) : {delta_max:.6f}")
#         print(f"Function Range (R)         : {r_val:.6f}")
#         print(f"Ratio (delta_max / R)      : {ratio:.6f}")
#         print("=" * 40)
#
#     if __name__ == "__main__":
#         main()
#
# if __name__ == "__main__":
#     main()

    """
    Numerical Analysis - Assignment Runner (Stage A Complete)
    Compares Neville and Cubic Spline interpolation, generates required graphs,
    tables, and displays polynomial coefficients.
    """


    def print_spline_coefficients(x_nodes, y_nodes):
        """
        Extracts and prints the a, b, c, d coefficients for each segment
        of the cubic spline, fulfilling requirement 3.
        """
        x, a = list(x_nodes), list(y_nodes)
        n = len(x)
        h = [x[i + 1] - x[i] for i in range(n - 1)]

        alpha, mu, z = [0] * n, [0] * n, [0] * n
        for i in range(1, n - 1):
            alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1])
            l = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
            mu[i] = h[i] / l
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l

        c = [0] * n
        for j in range(n - 2, -1, -1):
            c[j] = z[j] - mu[j] * c[j + 1]

        print("\n--- Spline Coefficients per Segment ---")
        print(f"{'Interval':<12} | {'a':<8} | {'b':<8} | {'c':<8} | {'d':<8}")
        print("-" * 55)
        for i in range(n - 1):
            b = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3
            d = (c[i + 1] - c[i]) / (3 * h[i])
            print(f"[{x[i]:.1f}, {x[i + 1]:.1f}]".ljust(12) + f" | {a[i]:>6.3f} | {b:>6.3f} | {c[i]:>6.3f} | {d:>6.3f}")


"""
Numerical Analysis - Assignment Runner (Stage A Complete - Standalone)
"""

import matplotlib.pyplot as plt
from typing import List, Tuple

Point = Tuple[float, float]


# ==========================================
# 1. Neville Interpolation Logic
# ==========================================
def neville_interpolation(points: List[Point], target_x: float) -> float:
    n = len(points)
    p = [point[1] for point in points]
    for j in range(1, n):
        for i in range(n - j):
            xi = points[i][0]
            xj = points[i + j][0]
            p[i] = ((target_x - xj) * p[i] - (target_x - xi) * p[i + 1]) / (xi - xj)
    return p[0]


# ==========================================
# 2. Cubic Spline Logic
# ==========================================
def cubic_spline(x_nodes, y_nodes, target_x):
    x, a = list(x_nodes), list(y_nodes)
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]

    alpha, mu, z = [0] * n, [0] * n, [0] * n
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1])
        l = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l

    c = [0] * n
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]

    i = 0
    while i < n - 2 and target_x > x[i + 1]:
        i += 1

    b = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3
    d = (c[i + 1] - c[i]) / (3 * h[i])
    dx = target_x - x[i]
    return a[i] + b * dx + c[i] * (dx ** 2) + d * (dx ** 3)


def print_spline_coefficients(x_nodes, y_nodes):
    """Calculates and prints the intermediate spline coefficients (Requirement 3)."""
    x, a = list(x_nodes), list(y_nodes)
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]

    alpha, mu, z = [0] * n, [0] * n, [0] * n
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1])
        l = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l

    c = [0] * n
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]

    print("\n--- Spline Coefficients per Segment ---")
    print(f"{'Interval':<12} | {'a':<8} | {'b':<8} | {'c':<8} | {'d':<8}")
    print("-" * 55)
    for i in range(n - 1):
        b = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3
        d = (c[i + 1] - c[i]) / (3 * h[i])
        print(f"[{x[i]:.1f}, {x[i + 1]:.1f}]".ljust(12) + f" | {a[i]:>6.3f} | {b:>6.3f} | {c[i]:>6.3f} | {d:>6.3f}")


# ==========================================
# 3. Main Assignment Runner
# ==========================================
def main():
    # Diagnostic print to ensure the script is running
    print("\n>>> Program Started. Calculating 200 points... <<<\n")

    # Original Data Points
    points = [
        (0.0, 0.500), (0.3, 0.779), (0.7, -0.054), (1.0, -0.591),
        (1.5, -0.368), (2.0, -0.715), (2.5, -0.056), (3.0, 0.976)
    ]
    x_nodes = [p[0] for p in points]
    y_nodes = [p[1] for p in points]

    # Generate 200 uniform points
    num_points = 200
    start_x, end_x = 0.0, 3.0
    test_points = [start_x + i * (end_x - start_x) / (num_points - 1) for i in range(num_points)]

    y_neville_list = []
    y_spline_list = []
    differences = []

    # Evaluate both methods over the 200 points
    for x in test_points:
        y_nev = neville_interpolation(points, x)
        y_spl = cubic_spline(x_nodes, y_nodes, x)

        y_neville_list.append(y_nev)
        y_spline_list.append(y_spl)
        differences.append(abs(y_nev - y_spl))

    # Calculate Metrics
    delta_max = max(differences)
    max_diff_index = differences.index(delta_max)
    x_delta_max = test_points[max_diff_index]

    r_val = max(y_nodes) - min(y_nodes)
    ratio = delta_max / r_val

    # Print Stage A Status
    print("=" * 55)
    print(" STAGE A: Validating the Continuous Function")
    print("=" * 55)
    print(f"Max Difference (delta_max) : {delta_max:.6f}")
    print(f"Function Range (R)         : {r_val:.6f}")
    print(f"Ratio (delta_max / R)      : {ratio:.6f}")

    print("\nStatus Check:")
    if ratio < 0.01:
        print(">> delta_max / R < 0.01 : Full agreement. Proceed to Stage B.")
    elif 0.01 <= ratio < 0.05:
        print(">> 0.01 <= delta_max / R < 0.05 : Partial agreement. Analyze and pass.")
    else:
        print(">> delta_max / R >= 0.05 : Discrepancy found. Stop, analyze, and fix.")

    # Requirement 3: Spline Coefficients & Neville at x=1.5
    print_spline_coefficients(x_nodes, y_nodes)
    val_1_5 = neville_interpolation(points, 1.5)
    print(f"\nExample: Value of Neville polynomial at x=1.5 is {val_1_5:.5f}")

    # Requirement 4: Table A
    table_x_vals = [0.5, 1.0, 1.5, 2.0, 2.5]
    print("\n--- Table A: Specific Point Comparisons ---")
    print(f"{'x':<5} | {'Neville (Y)':<15} | {'Spline (Y)':<15} | {'Difference (Delta)':<15}")
    print("-" * 55)
    for tx in table_x_vals:
        yn = neville_interpolation(points, tx)
        ys = cubic_spline(x_nodes, y_nodes, tx)
        print(f"{tx:<5.1f} | {yn:<15.6f} | {ys:<15.6f} | {abs(yn - ys):<15.6f}")
    print("=" * 55)

    # --- Plotting Requirement 1: Method Comparison ---
    plt.figure(figsize=(10, 5))
    plt.plot(test_points, y_neville_list, label='Neville Interpolation', color='blue', alpha=0.7)
    plt.plot(test_points, y_spline_list, label='Cubic Spline', color='red', linestyle='dashed')
    plt.scatter(x_nodes, y_nodes, color='black', zorder=5, label='Sampled Points')
    plt.title('Graph 1: Neville vs Cubic Spline Interpolation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()

    # --- Plotting Requirement 2: Delta Graph ---
    plt.figure(figsize=(10, 5))
    plt.plot(test_points, differences, label='delta(x) = |Neville - Spline|', color='purple')
    plt.scatter([x_delta_max], [delta_max], color='red', zorder=5)
    plt.annotate(f'delta_max\nx={x_delta_max:.2f}',
                 xy=(x_delta_max, delta_max), xytext=(x_delta_max, delta_max + 0.1),
                 arrowprops=dict(facecolor='black', shrink=0.05),
                 horizontalalignment='center')
    plt.title('Graph 2: Absolute Difference (Delta) Between Methods')
    plt.xlabel('x')
    plt.ylabel('Difference')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()