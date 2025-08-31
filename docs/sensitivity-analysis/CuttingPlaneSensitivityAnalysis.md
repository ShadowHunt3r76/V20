# CuttingPlaneSensitivityAnalysis

## Overview
`CuttingPlaneSensitivityAnalysis` provides sensitivity analysis for solutions obtained using the Cutting Plane algorithm, which is particularly useful for solving integer programming problems. This analysis helps understand how changes in the problem parameters affect the optimal integer solution and the cuts generated during the solution process.

## Key Features
- **LP-based Sensitivity Analysis**: Performs sensitivity analysis on the final LP relaxation
- **Cut-Specific Analysis**: Analyzes the impact of individual cuts on the solution
- **Objective Coefficient Ranges**: Determines allowable ranges for objective function coefficients
- **Constraint RHS Ranges**: Identifies valid ranges for constraint right-hand sides
- **Dual Values**: Provides shadow prices for both original constraints and cuts
- **Reduced Costs**: Calculates reduced costs for variables

## Usage

### Initialization
```csharp
var analysis = new CuttingPlaneSensitivityAnalysis(
    cuttingPlaneSolution,  // CuttingPlaneSolution - The solution to analyze
    originalModel         // ParsedLinearProgrammingModel - The original model
);
```

### Performing Analysis
```csharp
// Perform complete sensitivity analysis and print results
analysis.PerformAnalysis();

// Get specific sensitivity information
var objRange = analysis.GetObjectiveCoefficientRange(variableIndex);
var rhsRange = analysis.GetRightHandSideRange(constraintIndex);
var dualValues = analysis.GetDualValues();
var reducedCosts = analysis.GetReducedCosts();
bool isOptimal = analysis.IsOptimal();
```

## Example

### Problem Setup
Consider the following integer programming problem:
```
Maximize Z = 3x1 + 2x2
Subject to:
  2x1 + x2 ≤ 10
  x1 + x2 ≤ 6
  x1, x2 ≥ 0 and integer
```

### Performing Sensitivity Analysis
```csharp
// After solving the problem using Cutting Plane, we have:
var solution = new CuttingPlaneSolution
{
    BestObjectiveValue = 12.0,  // x1=2, x2=3
    BestSolution = new[] { 2.0, 3.0 },
    OverallStatus = SolutionStatus.Optimal,
    Iterations = 5,
    CutsAdded = 2,
    FinalLPSolution = new PrimalSimplexSolution { /* ... */ },
    Model = model,
    CutHistory = new[]
    {
        new CutInfo { Iteration = 1, CutCoefficients = new[] { -0.5, 1.0 }, Rhs = 1.5 },
        new CutInfo { Iteration = 3, CutCoefficients = new[] { 1.0, 0.0 }, Rhs = 2.0 }
    }
};

// Create and perform sensitivity analysis
var analysis = new CuttingPlaneSensitivityAnalysis(solution, model);
analysis.PerformAnalysis();
```

### Sample Output
```
==========================================================
CUTTING PLANE SENSITIVITY ANALYSIS
==========================================================

SOLUTION INFORMATION
----------------------------------------
Status: Optimal
Objective Value: 12.000000

Solution Vector:
x1: 2.000000 (Integer)
x2: 3.000000 (Integer)

Iterations: 5
Cuts Added: 2

LP-BASED SENSITIVITY ANALYSIS
----------------------------------------
[Standard LP sensitivity analysis output from the final LP relaxation...]

CUT-SPECIFIC SENSITIVITY ANALYSIS
----------------------------------------
Cut 1 (Iteration 1): -0.5x1 + x2 ≤ 1.5
  - Impact on Objective: -1.00
  - Dual Value: 0.50

Cut 2 (Iteration 3): x1 ≤ 2
  - Impact on Objective: -3.00
  - Dual Value: 1.00

CONSTRAINT SENSITIVITY
----------------------------------------
Constraint  Type     Slack    Dual Value   Status
1           ≤        1.00     0.00         Non-binding
2           ≤        1.00     0.00         Non-binding
Cut 1       ≤        0.50     0.50         Binding
Cut 2       ≤        0.00     1.00         Binding
```

## Interpretation of Results

### Cut-Specific Analysis
- **Cut 1**: The first cut added at iteration 1 removed part of the feasible region, reducing the objective by 1.00
- **Cut 2**: The second cut at iteration 3 was crucial, reducing the objective by 3.00 and making the solution integer-feasible

### Constraint Sensitivity
- The original constraints are non-binding (slack > 0)
- Both cuts are binding (slack = 0) with positive dual values, indicating they are actively constraining the solution
- The second cut has a higher dual value (1.00 vs 0.50), indicating it has a stronger impact on the objective

## Important Notes
1. **LP Relaxation**: The LP-based analysis is performed on the final LP relaxation, which includes all cuts
2. **Cut Impact**: The impact of cuts is measured by how much they reduce the objective value of the LP relaxation
3. **Dual Values**: For cuts, dual values indicate how much the objective would improve if the cut were relaxed by one unit
4. **Integer Solutions**: The analysis provides insights into the LP relaxation; the actual integer solution may behave differently

## See Also
- [SimplexSensitivityAnalysis](./SimplexSensitivityAnalysis.md)
- [RevisedSimplexSensitivityAnalysis](./RevisedSimplexSensitivityAnalysis.md)
- [BranchAndBoundSensitivityAnalysis](./BranchAndBoundSensitivityAnalysis.md)
- [BaseSensitivityAnalysis](./BaseSensitivityAnalysis.md)
