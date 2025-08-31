using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;
using static LinearProgramming.Algorithms.SensitivityAnalysis.SensitivityOutputFormatter;

namespace LinearProgramming.Algorithms.BranchAndBound
{
    /// <summary>
    /// Provides sensitivity analysis for Branch and Bound algorithm solutions
    /// </summary>
    public class BranchAndBoundSensitivityAnalysis : BaseSensitivityAnalysis
    {
        private readonly BranchAndBoundSolution _solution;
        private readonly ParsedLinearProgrammingModel _originalModel;
        private readonly SensitivityAnalysis _lpSensitivity;

        /// <summary>
        /// Initializes a new instance of the BranchAndBoundSensitivityAnalysis class
        /// </summary>
        /// <param name="solution">The Branch and Bound solution to analyze</param>
        /// <param name="originalModel">The original parsed model</param>
        public BranchAndBoundSensitivityAnalysis(BranchAndBoundSolution solution, ParsedLinearProgrammingModel originalModel)
            : base(
                originalModel?.Variables?.Select(v => v.Name).ToArray() ?? Array.Empty<string>(),
                originalModel?.ConstraintTypes ?? Array.Empty<ConstraintType>(),
                originalModel?.Variables?.Count ?? 0,
                originalModel?.ConstraintTypes?.Length ?? 0)
        {
            _solution = solution ?? throw new ArgumentNullException(nameof(solution));
            _originalModel = originalModel ?? throw new ArgumentNullException(nameof(originalModel));
            
            // Create LP sensitivity analysis for the best node's solution
            if (solution.BestCandidate?.Solution?.OptimalTable != null)
            {
                var basis = FindBasis(solution.BestCandidate.Solution.OptimalTable, 
                                   solution.CanonicalForm.ObjectiveCoefficients.Length);
                
                _lpSensitivity = new SensitivityAnalysis(
                    solution.BestCandidate.Solution.OptimalTable,
                    basis,
                    solution.BestCandidate.Model.Variables,
                    solution.BestCandidate.Model.ConstraintTypes,
                    solution.BestCandidate.Model.ObjectiveCoefficients.Length,
                    solution.BestCandidate.Model.ConstraintTypes.Length);
            }
        }

        /// <summary>
        /// Performs complete sensitivity analysis on the Branch and Bound solution
        /// </summary>
        public void PerformAnalysis()
        {
            if (_solution.BestCandidate == null)
            {
                Console.WriteLine("No optimal solution found for sensitivity analysis.");
                return;
            }

            PrintSectionHeader("BRANCH AND BOUND SENSITIVITY ANALYSIS");

            try
            {
                // 1. Basic solution information
                DisplaySolutionInfo();

                // 2. LP-based sensitivity analysis (from the best node)
                if (_lpSensitivity != null)
                {
                    PrintSubsectionHeader("LP-BASED SENSITIVITY ANALYSIS (FROM BEST NODE)");
                    PerformLPSensitivityAnalysis();
                }

                // 3. Integer-specific sensitivity analysis
                PrintSubsectionHeader("INTEGER-SPECIFIC SENSITIVITY ANALYSIS");
                PerformIntegerSensitivityAnalysis();
                
                PrintSectionHeader("END OF SENSITIVITY ANALYSIS");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nError during sensitivity analysis: {ex.Message}");
                if (ex.InnerException != null)
                {
                    Console.WriteLine($"Inner exception: {ex.InnerException.Message}");
                }
                Console.WriteLine("\n" + new string('=', 80));
            }
        }

        private void DisplaySolutionInfo()
        {
            PrintSubsectionHeader("SOLUTION INFORMATION");
            
            Console.WriteLine(FormatKeyValue("Status", _solution.OverallStatus));
            Console.WriteLine(FormatKeyValue("Objective Value", $"{_solution.BestObjectiveValue:F6}"));
            Console.WriteLine("\nSolution Vector:");
            
            for (int i = 0; i < _solution.BestSolution.Length; i++)
            {
                string varType = _originalModel.Variables[i].Type == VariableType.Binary ? " (Binary)" : 
                               _originalModel.Variables[i].Type == VariableType.Integer ? " (Integer)" : "";
                Console.WriteLine($"x{i + 1}{varType}: {_solution.BestSolution[i]:F6}");
            }
            
            Console.WriteLine($"\nTotal Nodes Explored: {_solution.TotalNodesExplored}");
            Console.WriteLine($"Nodes Fathomed: {_solution.TotalNodesFathomed}");
            Console.WriteLine($"Branching Depth: {_solution.BestCandidate?.Depth ?? 0}");
        }

        private void PerformLPSensitivityAnalysis()
        {
            try
            {
                // 1. Non-basic variable ranges
                Console.WriteLine("\nRANGES FOR NON-BASIC VARIABLES (AT OPTIMAL NODE):");
                for (int i = 0; i < _originalModel.Variables.Length; i++)
                {
                    if (!_solution.BestCandidate.Solution.BasisIndices.Contains(i))
                    {
                        try
                        {
                            var range = _lpSensitivity.GetNonBasicVariableRange(i);
                            Console.WriteLine($"Variable x{i + 1}: [{range.lowerBound:F4}, {range.upperBound:F4}]");
                        }
                        catch (Exception ex)
                        {
                            Console.WriteLine($"Error analyzing variable x{i + 1}: {ex.Message}");
                        }
                    }
                }

                // 2. Basic variable ranges
                Console.WriteLine("\nRANGES FOR BASIC VARIABLES (AT OPTIMAL NODE):");
                foreach (int varIndex in _solution.BestCandidate.Solution.BasisIndices)
                {
                    if (varIndex < _originalModel.Variables.Length) // Only original variables, not slacks
                    {
                        try
                        {
                            var range = _lpSensitivity.GetBasicVariableRange(varIndex);
                            Console.WriteLine($"Variable x{varIndex + 1}: [{range.lowerBound:F4}, {range.upperBound:F4}]");
                        }
                        catch (Exception ex)
                        {
                            Console.WriteLine($"Error analyzing basic variable x{varIndex + 1}: {ex.Message}");
                        }
                    }
                }

                // 3. Constraint RHS ranges
                Console.WriteLine("\nCONSTRAINT RHS RANGES (AT OPTIMAL NODE):");
                for (int i = 0; i < _solution.CanonicalForm.ConstraintTypes.Length; i++)
                {
                    try
                    {
                        var range = _lpSensitivity.GetRHSRange(i);
                        Console.WriteLine($"Constraint {i + 1} ({_solution.CanonicalForm.ConstraintTypes[i]}): [{range.lowerBound:F4}, {range.upperBound:F4}]");
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine($"Error analyzing constraint {i + 1}: {ex.Message}");
                    }
                }

                // 4. Shadow prices
                Console.WriteLine("\nSHADOW PRICES (DUAL VARIABLES):");
                for (int i = 0; i < _solution.CanonicalForm.ConstraintTypes.Length; i++)
                {
                    try
                    {
                        double shadowPrice = _lpSensitivity.GetShadowPrice(i);
                        Console.WriteLine($"Constraint {i + 1}: {shadowPrice:F4}");
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine($"Error getting shadow price for constraint {i + 1}: {ex.Message}");
                    }
                }

                // 5. Duality analysis
                try
                {
                    var dualityResult = _lpSensitivity.AnalyzeDuality();
                    Console.WriteLine("\nDUALITY ANALYSIS:");
                    Console.WriteLine($"Primal Objective: {dualityResult.PrimalObjective:F6}");
                    Console.WriteLine($"Dual Objective: {dualityResult.DualObjective:F6}");
                    Console.WriteLine($"Duality Gap: {dualityResult.DualityGap:E6}");
                    Console.WriteLine($"Strong Duality Holds: {dualityResult.HasStrongDuality}");
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"Error performing duality analysis: {ex.Message}");
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error performing LP-based sensitivity analysis: {ex.Message}");
            }
        }

        private void PerformIntegerSensitivityAnalysis()
        {
            try
            {
                // 1. Integer variable sensitivity
                Console.WriteLine("\nINTEGER VARIABLE SENSITIVITY:");
                for (int i = 0; i < _originalModel.Variables.Length; i++)
                {
                    if (_originalModel.Variables[i].Type == VariableType.Integer || 
                        _originalModel.Variables[i].Type == VariableType.Binary)
                    {
                        double value = _solution.BestSolution[i];
                        double floor = Math.Floor(value);
                        double ceil = Math.Ceiling(value);
                        
                        Console.WriteLine($"x{i + 1}: Value = {value:F0} (Integer)");
                        
                        // Check if we can increase/decrease by 1 and stay integer feasible
                        bool canIncrease = true;
                        bool canDecrease = true;
                        
                        // For binary variables, we can't go below 0 or above 1
                        if (_originalModel.Variables[i].Type == VariableType.Binary)
                        {
                            canIncrease = value < 1.0 - EPSILON;
                            canDecrease = value > EPSILON;
                            
                            Console.WriteLine("  Binary variable - " + 
                                $"Can {(canIncrease ? "" : "not ")}increase to 1, " +
                                $"Can {(canDecrease ? "" : "not ")}decrease to 0");
                        }
                        else
                        {
                            Console.WriteLine("  Integer variable - " + 
                                $"Can increase to {value + 1}, " +
                                $"Can decrease to {value - 1}");
                        }
                    }
                }

                // 2. Constraint tightness analysis
                Console.WriteLine("\nCONSTRAINT TIGHTNESS ANALYSIS:");
                var activeConstraints = new List<int>();
                var slackValues = new Dictionary<int, double>();
                
                for (int i = 0; i < _solution.CanonicalForm.ConstraintTypes.Length; i++)
                {
                    double lhs = 0;
                    for (int j = 0; j < _originalModel.Variables.Length; j++)
                    {
                        lhs += _solution.CanonicalForm.CoefficientMatrix[i][j] * _solution.BestSolution[j];
                    }
                    
                    double rhs = _solution.CanonicalForm.RHSVector[i];
                    double slack = 0;
                    
                    switch (_solution.CanonicalForm.ConstraintTypes[i])
                    {
                        case ConstraintType.LessThanOrEqual:
                            slack = rhs - lhs;
                            break;
                        case ConstraintType.GreaterThanOrEqual:
                            slack = lhs - rhs;
                            break;
                        case ConstraintType.Equal:
                            slack = Math.Abs(lhs - rhs);
                            break;
                    }
                    
                    slackValues[i] = slack;
                    
                    if (slack < EPSILON)
                    {
                        activeConstraints.Add(i);
                        Console.WriteLine($"Constraint {i + 1} is ACTIVE (slack = {slack:E6})");
                    }
                    else
                    {
                        Console.WriteLine($"Constraint {i + 1} has slack of {slack:F6}");
                    }
                }

                // 3. Branching variable impact
                if (_solution.BestCandidate != null && _solution.BestCandidate.BranchingConstraint != null)
                {
                    Console.WriteLine("\nBRANCHING IMPACT ANALYSIS:");
                    Console.WriteLine($"Last branching constraint: {_solution.BestCandidate.BranchingConstraint}");
                    
                    // Find siblings to compare impact of branching
                    var parentNode = _solution.AllNodes.FirstOrDefault(n => n.Id == _solution.BestCandidate.ParentId);
                    if (parentNode != null)
                    {
                        var siblings = _solution.AllNodes
                            .Where(n => n.ParentId == parentNode.Id && n.Id != _solution.BestCandidate.Id)
                            .ToList();
                        
                        if (siblings.Any())
                        {
                            Console.WriteLine("Alternative branches had the following outcomes:");
                            foreach (var sibling in siblings)
                            {
                                string status = sibling.IsFathomed ? 
                                    $"FATHOMED ({sibling.FathomReason})" : 
                                    "EXPLORED";
                                
                                Console.WriteLine($"- {sibling.BranchingConstraint}: {status}");
                                if (sibling.Solution?.Status == "Optimal")
                                {
                                    Console.WriteLine($"  Objective: {sibling.Solution.ObjectiveValue:F6}");
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error performing integer-specific sensitivity analysis: {ex.Message}");
            }
        }

        // Helper method to find basis variables from the final tableau with improved numerical stability
        private int[] FindBasis(double[,] tableau, int numOriginalVars)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            var basis = new List<int>();
            
            // First, find all columns that are part of the basis
            for (int j = 0; j < numOriginalVars; j++)
            {
                int basicRow = -1;
                bool isBasic = true;
                double maxValue = 0;

                // Check if this column is a unit vector (with some tolerance)
                for (int i = 1; i < rows; i++) // Skip objective row
                {
                    double absVal = Math.Abs(tableau[i, j]);
                    
                    // Check if this could be the '1' in a unit vector
                    if (absVal.AlmostEquals(1.0, EPSILON * 10)) // Slightly more tolerance for 1.0
                    {
                        if (basicRow == -1) 
                        {
                            basicRow = i - 1; // Convert to 0-based constraint index
                            maxValue = absVal;
                        }
                        else 
                        { 
                            // Found another non-zero in this column - not a unit vector
                            isBasic = false; 
                            break; 
                        }
                    }
                    else if (absVal > EPSILON)
                    {
                        // Non-zero value found where we expect zero
                        isBasic = false;
                        break;
                    }
                }

                // Additional check: make sure the row where we found the '1' doesn't have other non-zeros
                if (isBasic && basicRow != -1)
                {
                    for (int k = 0; k < cols; k++)
                    {
                        if (k != j && !tableau[basicRow + 1, k].IsZero(EPSILON))
                        {
                            isBasic = false;
                            break;
                        }
                    }
                }

                if (isBasic && basicRow != -1)
                {
                    basis.Add(j);
                }
            }
            
            // Verify we found exactly the right number of basic variables
            if (basis.Count != rows - 1) // rows-1 because we skip the objective row
            {
                Console.WriteLine($"⚠ Warning: Expected {rows-1} basis variables but found {basis.Count}. " +
                               "Results may be inaccurate due to numerical instability.");
            }
            
            return basis.ToArray();
        }
    }
}
