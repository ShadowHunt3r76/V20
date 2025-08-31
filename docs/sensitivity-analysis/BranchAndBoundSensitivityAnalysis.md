# BranchAndBoundSensitivityAnalysis

## Overview
`BranchAndBoundSensitivityAnalysis` provides sensitivity analysis for solutions obtained using the Branch and Bound algorithm, which is particularly useful for integer programming problems. This analysis helps understand how changes in the problem parameters affect the optimal integer solution.

## Key Features
- **LP-based Sensitivity Analysis**: Performs sensitivity analysis on the LP relaxation at the best node
- **Integer-Specific Analysis**: Provides insights specific to integer solutions
- **Objective Coefficient Ranges**: Determines allowable ranges for objective function coefficients
- **Constraint RHS Ranges**: Identifies valid ranges for constraint right-hand sides
- **Dual Values**: Provides shadow prices from the LP relaxation
- **Reduced Costs**: Calculates reduced costs for variables
- **Text-based Visualizations**: Includes rich text-based visualizations of sensitivity analysis results

## Usage

### Initialization
```csharp
var analysis = new BranchAndBoundSensitivityAnalysis(
    branchAndBoundSolution,  // BranchAndBoundSolution - The solution to analyze
    originalModel           // ParsedLinearProgrammingModel - The original model
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
Consider the following 0-1 Knapsack problem:
```
Maximize Z = 10x1 + 20x2 + 15x3
Subject to:
  5x1 + 10x2 + 7x3 ≤ 12
  x1, x2, x3 ∈ {0, 1}
```

### Performing Sensitivity Analysis
```csharp
// After solving the problem using Branch and Bound, we have:
var solution = new BranchAndBoundSolution
{
    BestObjectiveValue = 25.0,  // x1=1, x2=0, x3=1
    BestSolution = new[] { 1.0, 0.0, 1.0 },
    OverallStatus = SolutionStatus.Optimal,
    // ... other solution properties ...
};

var model = new ParsedLinearProgrammingModel
{
    Objective = new Objective { Type = OptimizationType.Maximize },
    Variables = new[]
    {
        new Variable { Name = "x1", Type = VariableType.Binary },
        new Variable { Name = "x2", Type = VariableType.Binary },
        new Variable { Name = "x3", Type = VariableType.Binary }
    },
    ObjectiveCoefficients = new[] { 10.0, 20.0, 15.0 },
    ConstraintCoefficients = new[] { new[] { 5.0, 10.0, 7.0 } },
    ConstraintConstants = new[] { 12.0 },
    ConstraintTypes = new[] { ConstraintType.LessThanOrEqual }
};

// Create and perform sensitivity analysis
var analysis = new BranchAndBoundSensitivityAnalysis(solution, model);
analysis.PerformAnalysis();
```

### Sample Output with Visualizations
```
==========================================================
BRANCH AND BOUND SENSITIVITY ANALYSIS
==========================================================

SOLUTION INFORMATION
----------------------------------------
Status: Optimal
Objective Value: 25.000000

Solution Vector:
x1: 1.000000 (Binary)
x2: 0.000000 (Binary)
x3: 1.000000 (Binary)

Total Nodes Explored: 5

LP-BASED SENSITIVITY ANALYSIS (FROM BEST NODE)
----------------------------------------
OBJECTIVE COEFFICIENT RANGES:
x1: [8.000000, 15.000000] (Current: 10.000000)
x2: [0.000000, 25.000000] (Current: 20.000000)
x3: [10.000000, 20.000000] (Current: 15.000000)

VISUALIZATION
----------------------------------------
OBJECTIVE COEFFICIENT RANGES
x1   [#############|------------------] 8.00 ──●── 15.00 (10.00)
x2   [------------------|#############] 0.00 ────────●── 25.00 (20.00)
x3   [###########|-------------------] 10.00 ──●── 20.00 (15.00)

NON-ZERO REDUCED COSTS
x2: ████████████████████████████████████████████ 5.00

RIGHT-HAND SIDE RANGES
Constraint 1 (≤) [10.00 ────●─── 15.00] (12.00)

SHADOW PRICES (DUAL VARIABLES)
Constraint 1 (≤): ██████████████████████████████████████████ 2.00

INTEGER-SPECIFIC SENSITIVITY ANALYSIS
----------------------------------------
VARIABLE  VALUE  TYPE     IMPACT ON OBJECTIVE
x1        1      Binary   Decrease if excluded: -10.00
x2        0      Binary   Increase if included: +20.00 (Infeasible)
x3        1      Binary   Decrease if excluded: -15.00

CONSTRAINT SENSITIVITY
----------------------------------------
CONSTRAINT  TYPE           SLACK    DUAL VALUE
1           ≤ (Knapsack)   0.00     2.00
```

## Interpretation of Results

### Visualizations Explained
1. **Objective Coefficient Ranges**:
   - Bar charts show the allowable range for each objective coefficient
   - The `●` marker indicates the current coefficient value
   - Example: `x1` can range from 8.00 to 15.00 (current: 10.00)

2. **Reduced Costs**:
   - Histogram shows non-zero reduced costs for non-basic variables
   - Indicates how much the objective would change if the variable entered the basis

3. **Right-Hand Side Ranges**:
   - Shows the valid range for each constraint's RHS
   - The `●` marker indicates the current RHS value

4. **Shadow Prices**:
   - Histogram shows the dual value for each binding constraint
   - Indicates the rate of change of the objective per unit change in the RHS

### Integer-Specific Analysis
- **x1 (Value: 1)**: Excluding x1 would decrease the objective by 10.00 (its coefficient)
- **x2 (Value: 0)**: Including x2 would increase the objective by 20.00, but is infeasible due to knapsack capacity
- **x3 (Value: 1)**: Excluding x3 would decrease the objective by 15.00

### Constraint Sensitivity
- The knapsack constraint is binding (slack = 0)
- The dual value of 2.00 indicates that increasing the knapsack capacity by 1 unit would improve the objective by 2.00 (in the LP relaxation)
- The RHS range shows that this rate holds for capacity between 10.00 and 15.00

## Important Notes
1. **LP Relaxation vs. Integer Solution**: The LP-based analysis is performed on the relaxation at the best node, which may differ from the integer solution.
2. **Discrete Nature**: For integer problems, the actual change in the optimal solution when parameters change may be different from the continuous case.
3. **Multiple Optimal Solutions**: The analysis corresponds to the specific optimal solution found; other optimal solutions may exist.
4. **Infeasibility**: Some parameter changes that would be valid in the LP relaxation might lead to integer infeasibility.

## See Also
- [SimplexSensitivityAnalysis](./SimplexSensitivityAnalysis.md)
- [RevisedSimplexSensitivityAnalysis](./RevisedSimplexSensitivityAnalysis.md)
- [BaseSensitivityAnalysis](./BaseSensitivityAnalysis.md)
