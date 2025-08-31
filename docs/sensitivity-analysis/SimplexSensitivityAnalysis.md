# SimplexSensitivityAnalysis

## Overview
`SimplexSensitivityAnalysis` is a class that performs sensitivity analysis on the optimal solution obtained from the Simplex algorithm. It helps understand how changes in the problem's coefficients affect the optimal solution, providing valuable insights for decision-making.

## Key Features
- **Objective Coefficient Ranges**: Determines the allowable range for each objective coefficient while maintaining the current optimal basis.
- **Right-Hand Side Ranges**: Identifies the range of values for each constraint's right-hand side that keeps the current basis optimal.
- **Dual Values**: Provides the dual values (shadow prices) for each constraint.
- **Reduced Costs**: Calculates the reduced costs for non-basic variables.
- **Optimality Verification**: Verifies if the current solution remains optimal under coefficient changes.

## Usage

### Initialization
```csharp
var sensitivityAnalysis = new SimplexSensitivityAnalysis(
    finalTableau,      // double[,] - Final simplex tableau
    basisIndices,      // int[] - Indices of basic variables
    variableNames,     // string[] - Names of all variables
    objectiveCoeffs,   // double[] - Original objective coefficients
    rhsCoefficients,   // double[] - Right-hand side coefficients
    constraintTypes    // ConstraintType[] - Types of constraints
);
```

### Performing Analysis
```csharp
// Perform complete sensitivity analysis and print results
sensitivityAnalysis.PerformAnalysis();

// Get specific sensitivity information
var objRange = sensitivityAnalysis.GetObjectiveCoefficientRange(variableIndex);
var rhsRange = sensitivityAnalysis.GetRightHandSideRange(constraintIndex);
var dualValues = sensitivityAnalysis.GetDualValues();
var reducedCosts = sensitivityAnalysis.GetReducedCosts();
bool isOptimal = sensitivityAnalysis.IsOptimal();
```

## Example

### Problem Setup
Consider the following linear programming problem:
```
Maximize Z = 3x + 2y
Subject to:
  2x + y ≤ 100
   x + y ≤ 80
   x ≤ 40
   x, y ≥ 0
```

### Performing Sensitivity Analysis
```csharp
// After solving the problem using Simplex, we have:
double[,] finalTableau = {
    // Z    x     y     s1    s2    s3    RHS
    {1.0,  0.0,  0.0,  1.0,  1.0,  0.0,  180.0},  // Z equation
    {0.0,  0.0,  1.0,  1.0, -1.0,  0.0,  20.0},   // y equation
    {0.0,  0.0,  0.0, -1.0,  2.0,  1.0,  20.0},   // s3 equation
    {0.0,  1.0,  0.0,  0.0,  1.0,  0.0,  40.0}    // x equation
};

int[] basisIndices = { -1, 1, 4, 0 };  // Z, y, s3, x
string[] variableNames = { "x", "y", "s1", "s2", "s3" };
double[] objectiveCoeffs = { 3.0, 2.0, 0.0, 0.0, 0.0 };
double[] rhsCoefficients = { 100.0, 80.0, 40.0 };
var constraintTypes = new[] { 
    ConstraintType.LessThanOrEqual, 
    ConstraintType.LessThanOrEqual, 
    ConstraintType.LessThanOrEqual 
};

// Create and perform sensitivity analysis
var analysis = new SimplexSensitivityAnalysis(
    finalTableau,
    basisIndices,
    variableNames,
    objectiveCoeffs,
    rhsCoefficients,
    constraintTypes
);

analysis.PerformAnalysis();
```

### Sample Output
```
==========================================================
SIMPLEX SENSITIVITY ANALYSIS
==========================================================

SOLUTION INFORMATION
----------------------------------------
Status: Optimal
Objective Value: 180.000000

VARIABLE INFORMATION
----------------------------------------
Name    Value    Status     Reduced Cost
x       40.00    Basic      0.00
y       20.00    Basic      0.00
s1       0.00    Non-Basic  1.00
s2       0.00    Non-Basic  1.00
s3      20.00    Basic      0.00

OBJECTIVE COEFFICIENT RANGES
----------------------------------------
Variable  Current Coef  Allowable Decrease  Allowable Increase
x         3.00          1.00                Infinity
y         2.00          1.00                1.00
s1        0.00          -Infinity           1.00
s2        0.00          -Infinity           1.00
s3        0.00          -1.00               Infinity

RIGHT-HAND SIDE RANGES
----------------------------------------
Constraint  Current RHS  Allowable Decrease  Allowable Increase
1           100.00       20.00                20.00
2            80.00       60.00                20.00
3            40.00       20.00                Infinity

DUAL VALUES
----------------------------------------
Constraint  Dual Value
1           1.00
2           1.00
3           0.00
```

## Interpretation of Results

### Objective Coefficient Ranges
- The coefficient of `x` in the objective function can decrease by 1.0 or increase without bound without changing the optimal basis.
- The coefficient of `y` can vary between 1.0 and 3.0 while maintaining optimality.
- The reduced costs show how much the objective would change if a non-basic variable were to enter the basis.

### Right-Hand Side Ranges
- The first constraint's RHS can vary between 80 and 120 without changing the optimal basis.
- The dual value of 1.0 for the first constraint indicates that increasing its RHS by 1 unit would improve the objective by 1 unit (within the allowable range).

## Notes
- The analysis assumes that only one coefficient changes at a time.
- The ranges are valid only within the specified limits; beyond these limits, the current basis may no longer be optimal.
- The analysis is based on the final simplex tableau, so the problem must be solved to optimality first.

## See Also
- [RevisedSimplexSensitivityAnalysis](./RevisedSimplexSensitivityAnalysis.md)
- [BaseSensitivityAnalysis](./BaseSensitivityAnalysis.md)
- [NumericalStabilityUtils](../utils/NumericalStabilityUtils.md)
