"""
Numerical Analysis - Interpolation Methods
Assignment: Polynomial Interpolation
Student: 2nd Year, Software Engineering
"""

import math


def polynomial_interpolation(table: list[tuple[float, float]], x: float) -> float:
    """
    Estimate the value of an unknown function at point x using polynomial
    interpolation (Lagrange basis form).

    How it works
    ------------
    Given n data points (x0,y0), (x1,y1), ..., (x_{n-1}, y_{n-1}), we build
    the unique polynomial P of degree (n-1) that passes through all of them.

    We write P as a weighted sum of n "basis polynomials" L_i(x):

        P(x) = y0*L0(x) + y1*L1(x) + ... + y_{n-1}*L_{n-1}(x)

    Each basis polynomial L_i(x) is designed so that:
        - L_i(x_i) = 1   (evaluates to 1 at its own x-node)
        - L_i(x_j) = 0   for every j ≠ i  (evaluates to 0 at all other nodes)

    This guarantees that P(x_i) = y_i for every data point — the polynomial
    passes exactly through all given points.

    L_i is built as a product over all nodes except x_i:

        L_i(x) = product of (x - x_j) / (x_i - x_j)  for all j ≠ i

    Args:
        table : List of (x, y) tuples representing the known data points.
                All x-values must be distinct.
                Example: [(1.0, 2.0), (2.0, 5.0), (3.0, 10.0)]
        x     : The x-value at which to estimate the function.

    Returns:
        The interpolated y-value P(x).

    Raises:
        ValueError: If fewer than 2 points are given, or if any x-values repeat.
    """

    # ── Validation ────────────────────────────────────────────────────────────
    if len(table) < 2:
        raise ValueError("At least 2 data points are required for interpolation.")

    # Unpack the table into separate x and y lists for easy indexing.
    x_points = [point[0] for point in table]
    y_points = [point[1] for point in table]

    # Check for duplicate x-values — they would cause a division by zero below
    # and also make the interpolation problem unsolvable (two y's for one x).
    if len(set(x_points)) != len(x_points):
        raise ValueError("All x-values in the table must be distinct.")

    # ── Shortcut: x is already a known point ─────────────────────────────────
    # No computation needed — we already have the exact answer.
    if x in x_points:
        return y_points[x_points.index(x)]

    # ── Core Lagrange interpolation ───────────────────────────────────────────
    n      = len(x_points)
    result = 0.0  # This will accumulate P(x) = sum of y_i * L_i(x)

    for i in range(n):
        # Build the i-th Lagrange basis polynomial evaluated at x.
        # Start with L_i(x) = 1 and multiply in each factor one by one.
        L_i = 1.0

        for j in range(n):
            if j == i:
                # Skip the i-th node itself — L_i has no (x - x_i) factor.
                continue

            # Multiply by the j-th factor: (x - x_j) / (x_i - x_j)
            # Numerator   (x - x_j)  : shifts the polynomial so it has a root at x_j
            # Denominator (x_i - x_j): scales it so L_i(x_i) = 1
            L_i *= (x - x_points[j]) / (x_points[i] - x_points[j])

        # Add y_i * L_i(x) to the running total.
        # Because L_i(x_k) = 0 for k≠i and L_i(x_i) = 1, each term contributes
        # exactly y_i when evaluated at its own node and 0 at all other nodes.
        result += y_points[i] * L_i

    return result

def linear_interpolation(table, x_target):
    """
    פונקציה לביצוע אינטרפולציה לינארית.
    מקבלת טבלה של נקודות (רשימה של טאפלים (x, y)) ונקודת יעד (x_target).
    מחזירה את ערך ה-Y המחושב.
    """
    # מיון הטבלה לפי ערכי ה-X בסדר עולה כדי להבטיח מציאת טווח נכון
    sorted_table = sorted(table, key=lambda point: point[0])
    
    # חיפוש שתי הנקודות התוחמות את נקודת היעד (x_target)
    for i in range(len(sorted_table) - 1):
        x1, y1 = sorted_table[i]
        x2, y2 = sorted_table[i + 1]
        
        # בדיקה האם ה-X המבוקש נמצא בין שתי הנקודות הנוכחיות
        if x1 <= x_target <= x2:
            # נוסחת אינטרפולציה לינארית: y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
            y_target = y1 + ((x_target - x1) * (y2 - y1) / (x2 - x1))
            return y_target

    # במקרה שהנקודה מחוץ לטווח הטבלה (אקסטרפולציה)
    # נשתמש בשתי הנקודות הראשונות או האחרונות בהתאמה
    if x_target < sorted_table[0][0]:
        x1, y1 = sorted_table[0]
        x2, y2 = sorted_table[1]
    else:
        x1, y1 = sorted_table[-2]
        x2, y2 = sorted_table[-1]
        
    y_target = y1 + ((x_target - x1) * (y2 - y1) / (x2 - x1))
    return y_target


if __name__ == "__main__":
    print("=" * 60)
    print("    Numerical Analysis - Interpolation Comparison")
    print("=" * 60)

    # 1. אינטרפולציה ליניארית (מבוצעת ראשונה)
    table_points = [(1, 2.5), (2, 4.0), (4, 7.5), (5, 9.0)]
    point_to_find = 3

    linear_result = linear_interpolation(table_points, point_to_find)

    print("\n[1] Linear Interpolation Result:")
    print(f"  Table        : {table_points}")
    print(f"  Query x      : {point_to_find}")
    print(f"  Result (Y)   : {linear_result:.6f}")

    # 2. אינטרפולציה פולינומיאלית על אותה הטבלה לשם השוואה
    poly_compare_result = polynomial_interpolation(table_points, point_to_find)
    print("\n[2] Polynomial Interpolation (Lagrange) on the same Table:")
    print(f"  Result (Y)   : {poly_compare_result:.6f}")

    # 3. דוגמאות נוספות של אינטרפולציה פולינומיאלית
    print("\n" + "-" * 60)
    print("    Additional Polynomial Interpolation Examples")
    print("-" * 60)

    # דוגמה א' — f(x) = x² + 1
    table_1 = [
        (1.0,  2.0),
        (2.0,  5.0),
        (3.0, 10.0),
        (4.0, 17.0),
    ]
    x1 = 2.5
    exact1 = x1**2 + 1
    result1 = polynomial_interpolation(table_1, x1)

    print(f"\nExample 1 — f(x) = x² + 1,  query x = {x1}")
    print(f"  Interpolated : {result1:.6f}")
    print(f"  Exact        : {exact1:.6f}")
    print(f"  Error        : {abs(result1 - exact1):.2e}")

    # דוגמה ב' — f(x) = sin(x)
    table_2 = [(xi, math.sin(xi)) for xi in [0.0, 0.5, 1.0, 1.5, 2.0]]
    x2 = 0.75
    exact2 = math.sin(x2)
    result2 = polynomial_interpolation(table_2, x2)

    print(f"\nExample 2 — f(x) = sin(x),  query x = {x2}")
    print(f"  Interpolated : {result2:.6f}")
    print(f"  Exact        : {exact2:.6f}")
    print(f"  Error        : {abs(result2 - exact2):.2e}")

    print("\n" + "=" * 60)