# RevisedSimplexSensitivityAnalysis

## Overview
`RevisedSimplexSensitivityAnalysis` extends the functionality of `BaseSensitivityAnalysis` to provide sensitivity analysis specifically for solutions obtained using the Revised Simplex method. This implementation is more memory-efficient for large-scale problems as it works with the basis inverse rather than the full tableau.

## Key Features
- **Efficient Analysis**: Uses the basis inverse for memory efficiency with large problems
- **Objective Coefficient Ranges**: Determines allowable ranges for objective function coefficients
- **Right-Hand Side Ranges**: Identifies valid ranges for constraint right-hand sides
- **Dual Values**: Provides shadow prices for constraints
- **Reduced Costs**: Calculates reduced costs for non-basic variables
- **Optimality Verification**: Checks if the current solution remains optimal

## Usage

### Initialization
```csharp
var sensitivityAnalysis = new RevisedSimplexSensitivityAnalysis(
    basisInverse,       // double[,] - Basis inverse matrix
    basicSolution,      // double[] - Basic solution values
    reducedCosts,       // double[] - Reduced costs for all variables
    basicVariables,     // int[] - Indices of basic variables
    variableNames,      // string[] - Names of all variables
    objectiveCoeffs,    // double[] - Original objective coefficients
    rhsCoefficients,    // double[] - Right-hand side coefficients
    constraintTypes     // ConstraintType[] - Types of constraints
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
// After solving the problem using Revised Simplex, we have:
double[,] basisInverse = {
    { 1.0,  0.0,  0.0 },
    { 0.0,  1.0,  0.0 },
    { 0.0,  0.0,  1.0 }
};

double[] basicSolution = { 180.0, 20.0, 20.0, 40.0 }; // Z, y, s3, x
double[] reducedCosts = { 0.0, 0.0, 1.0, 1.0, 0.0 };  // x, y, s1, s2, s3
int[] basicVariables = { -1, 1, 4, 0 };  // Z, y, s3, x
string[] variableNames = { "x", "y", "s1", "s2", "s3" };
double[] objectiveCoeffs = { 3.0, 2.0, 0.0, 0.0, 0.0 };
double[] rhsCoefficients = { 100.0, 80.0, 40.0 };
var constraintTypes = new[] { 
    ConstraintType.LessThanOrEqual, 
    ConstraintType.LessThanOrEqual, 
    ConstraintType.LessThanOrEqual 
};

// Create and perform sensitivity analysis
var analysis = new RevisedSimplexSensitivityAnalysis(
    basisInverse,
    basicSolution,
    reducedCosts,
    basicVariables,
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
REVISED SIMPLEX SENSITIVITY ANALYSIS
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

## Advantages Over Standard Simplex Sensitivity Analysis
- More memory efficient for large problems
- Better numerical stability due to working with the basis inverse
- Faster updates when only small changes are made to the problem

## Notes
- The analysis assumes that only one coefficient changes at a time.
- The ranges are valid only within the specified limits.
- The problem must be solved to optimality before performing sensitivity analysis.
- For maximization problems, non-basic variables must have non-negative reduced costs for optimality.
- For minimization problems, non-basic variables must have non-positive reduced costs for optimality.

## See Also
- [SimplexSensitivityAnalysis](./SimplexSensitivityAnalysis.md)
- [BaseSensitivityAnalysis](./BaseSensitivityAnalysis.md)
- [NumericalStabilityUtils](../utils/NumericalStabilityUtils.md)
