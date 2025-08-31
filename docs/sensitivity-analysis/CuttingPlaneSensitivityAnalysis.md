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
- **Text-based Visualizations**: Includes rich text-based visualizations of sensitivity analysis results

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

### Sample Output with Visualizations
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
OBJECTIVE COEFFICIENT RANGES:
x1: [2.000000, 4.000000] (Current: 3.000000)
x2: [1.500000, 3.000000] (Current: 2.000000)

VISUALIZATION
----------------------------------------
OBJECTIVE COEFFICIENT RANGES
x1   [#############|------------------] 2.00 ──●── 4.00 (3.00)
x2   [###########|-------------------] 1.50 ──●── 3.00 (2.00)

NON-ZERO REDUCED COSTS
[No non-basic variables with non-zero reduced costs]

RIGHT-HAND SIDE RANGES
Constraint 1 (≤) [8.00 ────●─── 12.00] (10.00)
Constraint 2 (≤) [4.00 ────●─── 8.00] (6.00)

SHADOW PRICES (DUAL VARIABLES)
Constraint 1 (≤): ██████████████████████████████████████████ 0.00
Constraint 2 (≤): ██████████████████████████████████████████ 0.00
Cut 1 (≤):        ██████████████████████████████████████████ 0.50
Cut 2 (≤):        ██████████████████████████████████████████ 1.00

CUT-SPECIFIC SENSITIVITY ANALYSIS
----------------------------------------
CUT 1 (Iteration 1): -0.5x1 + x2 ≤ 1.5
  - Impact on Objective: -1.00
  - Dual Value: 0.50
  - Slack: 0.50
  - Status: Binding

CUT 2 (Iteration 3): x1 ≤ 2
  - Impact on Objective: -3.00
  - Dual Value: 1.00
  - Slack: 0.00
  - Status: Tightly Binding

CONSTRAINT SENSITIVITY
----------------------------------------
CONSTRAINT  TYPE     SLACK    DUAL VALUE   STATUS
1           ≤        1.00     0.00         Non-binding
2           ≤        1.00     0.00         Non-binding
Cut 1       ≤        0.50     0.50         Binding
Cut 2       ≤        0.00     1.00         Tightly Binding
```

## Interpretation of Results

### Visualizations Explained
1. **Objective Coefficient Ranges**:
   - Bar charts show the allowable range for each objective coefficient
   - The `●` marker indicates the current coefficient value
   - Example: `x1` can range from 2.00 to 4.00 (current: 3.00)

2. **Right-Hand Side Ranges**:
   - Shows the valid range for each constraint's RHS
   - The `●` marker indicates the current RHS value
   - Example: Constraint 1's RHS can range from 8.00 to 12.00 (current: 10.00)

3. **Shadow Prices**:
   - Histogram shows the dual value for each constraint and cut
   - The length of the bar indicates the relative strength of the dual value
   - Cuts have positive dual values, showing their impact on the solution

### Cut-Specific Analysis
- **Cut 1 (Iteration 1)**:
  - Slightly binding (slack = 0.50)
  - Dual value of 0.50 indicates the objective would improve by 0.5 if the RHS increased by 1
  - Reduced the objective by 1.00 from the initial LP relaxation

- **Cut 2 (Iteration 3)**:
  - Tightly binding (slack = 0.00)
  - Higher dual value (1.00) shows it's more restrictive than Cut 1
  - Reduced the objective by 3.00, making the solution integer-feasible

### Constraint Sensitivity
- **Original Constraints**:
  - Both constraints are non-binding (slack = 1.00)
  - Dual values are 0.00, indicating they don't constrain the current solution

- **Cuts**:
  - Both cuts are binding (slack ≤ 0.50)
  - The second cut is more restrictive (higher dual value)
  - The dual values show the rate of objective improvement per unit relaxation of each cut

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
