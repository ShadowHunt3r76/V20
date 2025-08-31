# KnapsackSensitivityAnalysis

## Overview
`KnapsackSensitivityAnalysis` provides specialized sensitivity analysis for solutions obtained from the 0-1 Knapsack problem. It helps understand how changes in item values, weights, and knapsack capacity affect the optimal solution.

## Key Features
- **Item Value Sensitivity**: Analyzes how much item values can change without affecting the optimal solution
- **Weight Sensitivity**: Determines how changes in item weights impact the solution
- **Capacity Sensitivity**: Shows how the solution changes with different knapsack capacities
- **Critical Items**: Identifies items that are on the "borderline" of being included/excluded
- **Reduced Costs**: Provides reduced costs for items not in the optimal solution

## Usage

### Initialization
```csharp
var analysis = new KnapsackSensitivityAnalysis(
    knapsackSolution,  // KnapsackSolution - The solution to analyze
    originalModel     // ParsedLinearProgrammingModel - The original model
);
```

### Performing Analysis
```csharp
// Perform complete sensitivity analysis and print results
analysis.PerformAnalysis();

// Get specific sensitivity information
var valueRange = analysis.GetObjectiveCoefficientRange(itemIndex);
var capacityRange = analysis.GetRightHandSideRange(0);  // Only one constraint in knapsack
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
  5x1 + 10x2 + 7x3 ≤ 12  (knapsack capacity)
  x1, x2, x3 ∈ {0, 1}
```

### Performing Sensitivity Analysis
```csharp
// After solving the problem, we have:
var solution = new KnapsackSolution
{
    BestObjectiveValue = 25.0,  // x1=1, x2=0, x3=1 (value=25, weight=12)
    BestSolution = new[] { 1.0, 0.0, 1.0 },
    OverallStatus = SolutionStatus.Optimal,
    TotalItems = 3,
    KnapsackCapacity = 12,
    Items = new[]
    {
        new KnapsackItem { Value = 10, Weight = 5, IsSelected = true },
        new KnapsackItem { Value = 20, Weight = 10, IsSelected = false },
        new KnapsackItem { Value = 15, Weight = 7, IsSelected = true }
    }
};

// Create and perform sensitivity analysis
var analysis = new KnapsackSensitivityAnalysis(solution, model);
analysis.PerformAnalysis();
```

### Sample Output
```
==========================================================
KNAPSACK SENSITIVITY ANALYSIS
==========================================================

SOLUTION INFORMATION
----------------------------------------
Status: Optimal
Objective Value: 25.000000
Knapsack Capacity: 12.00
Total Weight: 12.00 (100.00% of capacity)

ITEM SENSITIVITY
----------------------------------------
Item  Value  Weight  Selected  Value Range        Weight Range
----  -----  ------  --------  ------------------  ------------------
x1    10.00   5.00    Yes      5.00 ≤ c₁ ≤ 15.00  4.00 ≤ w₁ ≤ 12.00
x2    20.00  10.00    No      15.00 ≤ c₂ < 20.00  7.00 < w₂ ≤ 15.00
x3    15.00   7.00    Yes     10.00 ≤ c₃ ≤ 20.00  5.00 ≤ w₃ ≤ 12.00

CAPACITY SENSITIVITY
----------------------------------------
Current Capacity: 12.00
Allowable Decrease: 2.00 (to 10.00)
Allowable Increase: 3.00 (to 15.00)

DUAL VALUES
----------------------------------------
Constraint  Dual Value  Status
Capacity     1.25        Binding

REDUCED COSTS
----------------------------------------
Item  Reduced Cost  Interpretation
----  ------------  ---------------------------------
x1     -5.00         Value can decrease by 5.00 before x1 is excluded
x2      5.00        Value must increase by 5.00 for x2 to be included
x3     -5.00        Value can decrease by 5.00 before x3 is excluded
```

## Interpretation of Results

### Item Sensitivity
- **x1 (Included)**: 
  - Value can decrease to 5.00 or increase to 15.00 without changing the solution
  - Weight can increase to 12.00 or decrease to 4.00 without changing the solution
- **x2 (Excluded)**:
  - Value needs to increase to 25.00 to be included
  - Weight needs to decrease to 7.00 to be included
- **x3 (Included)**:
  - Value can decrease to 10.00 or increase to 20.00 without changing the solution
  - Weight can increase to 12.00 or decrease to 5.00 without changing the solution

### Capacity Sensitivity
- The current capacity is 12.00
- The solution remains optimal if capacity decreases by up to 2.00 (to 10.00) or increases by up to 3.00 (to 15.00)

### Dual Value
- The dual value of 1.25 for the capacity constraint indicates that increasing the knapsack capacity by 1 unit would improve the objective by 1.25 units (within the allowable range)

### Reduced Costs
- Negative reduced costs for included items show how much their values can decrease before being excluded
- Positive reduced costs for excluded items show how much their values must increase to be included

## Important Notes
1. **Discrete Nature**: The analysis provides ranges within which the current solution remains optimal, but the actual change might lead to a different integer solution
2. **Single Parameter Changes**: The analysis assumes only one parameter changes at a time
3. **Ties**: When the solution is on the boundary of changing, the analysis indicates the critical values where the solution might change
4. **Infeasibility**: Some parameter changes might make the problem infeasible (e.g., reducing capacity below the weight of the lightest item)

## See Also
- [SimplexSensitivityAnalysis](./SimplexSensitivityAnalysis.md)
- [BranchAndBoundSensitivityAnalysis](./BranchAndBoundSensitivityAnalysis.md)
- [BaseSensitivityAnalysis](./BaseSensitivityAnalysis.md)
