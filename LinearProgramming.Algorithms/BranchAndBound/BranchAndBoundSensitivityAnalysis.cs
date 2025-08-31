using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using SensitivityAnalysisImpl = LinearProgramming.Algorithms.SensitivityAnalysis.SensitivityAnalysisImpl;
using LinearProgramming.Algorithms;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;
using static LinearProgramming.Algorithms.SensitivityAnalysis.SensitivityOutputFormatter;
using static LinearProgramming.Algorithms.PrimalSimplex.MatrixUtils;

namespace LinearProgramming.Algorithms.BranchAndBound
{
    /// <summary>
    /// Provides sensitivity analysis for Branch and Bound algorithm solutions
    /// </summary>
    public class BranchAndBoundSensitivityAnalysis : BaseSensitivityAnalysis, ISensitivityAnalysis
    {
        private readonly BranchAndBoundSolution _solution;
        private readonly ParsedLinearProgrammingModel _originalModel;
        private readonly SensitivityAnalysisImpl? _lpSensitivity;
        
        // Visualization data storage
        private Dictionary<string, (double min, double current, double max)> _objectiveRanges = new();
        private Dictionary<string, (double min, double current, double max)> _rhsRanges = new();
        private Dictionary<string, double> _shadowPrices = new();
        private Dictionary<string, double> _reducedCosts = new();

        /// <summary>
        /// Initializes a new instance of the BranchAndBoundSensitivityAnalysis class
        /// </summary>
        /// <param name="solution">The Branch and Bound solution to analyze</param>
        /// <param name="originalModel">The original parsed model</param>
        public BranchAndBoundSensitivityAnalysis(BranchAndBoundSolution solution, ParsedLinearProgrammingModel originalModel)
            : base(
                originalModel?.Variables?.Select(v => v.Name).ToArray() ?? Array.Empty<string>(),
                originalModel?.Constraints?.Select(c => c.Type).ToArray() ?? Array.Empty<ConstraintType>(),
                originalModel?.Variables?.Count ?? 0,
                originalModel?.Constraints?.Count ?? 0)
        {
            _solution = solution ?? throw new ArgumentNullException(nameof(solution));
            _originalModel = originalModel ?? throw new ArgumentNullException(nameof(originalModel));
            
            // Create LP sensitivity analysis for the best node's solution
            if (solution.BestCandidate?.Solution?.OptimalTable != null)
            {
                var basis = FindBasis(solution.BestCandidate.Solution.OptimalTable, 
                                   _originalModel.Variables.Count); // Use original variable count
                
                var canonicalModel = solution.BestCandidate.Model;
                var optimalTable = solution.BestCandidate.Solution.OptimalTable;
                
                // Create variables array with correct types and coefficients
                var variables = new List<LinearProgramming.Algorithms.Variable>();
                for (int i = 0; i < canonicalModel.VariableTypes.Length; i++)
                {
variables.Add(new LinearProgramming.Algorithms.Variable 
                    { 
                        Name = i < _originalModel.Variables.Count ? _originalModel.Variables[i].Name : $"s{i}",
                        Coefficient = i < canonicalModel.ObjectiveCoefficients.Length ? 
                            canonicalModel.ObjectiveCoefficients[i] : 0,
                        LowerBound = 0,
                        UpperBound = double.PositiveInfinity,
                        IsInteger = i < _originalModel.Variables.Count ? 
                            (_originalModel.Variables[i].Type == VariableType.Integer || 
                             _originalModel.Variables[i].Type == VariableType.Binary) : false
                    });
                }
                
                // Initialize the LP sensitivity analysis implementation with the correct variable types
                _lpSensitivity = new SensitivityAnalysisImpl(
                    optimalTable,
                    basis,
                    variables.ToArray(),
                    canonicalModel.ConstraintTypes,
                    _originalModel.Variables.Count, // Number of original variables
                    canonicalModel.ConstraintTypes.Length);
            }
        }

        /// <summary>
        /// Performs complete sensitivity analysis on the Branch and Bound solution
        /// </summary>
        public override void PerformAnalysis()
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
                
                // 4. Generate and display visualization
                if (_lpSensitivity != null)
                {
                    PrintSubsectionHeader("VISUALIZATION");
                    Console.WriteLine(GenerateVisualization());
                }
                
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
                if (_solution.BestCandidate?.Solution is ILinearProgramSolution lpSolution && 
                    _lpSensitivity != null &&
                    _originalModel?.Variables != null)
                {
                    // Clear previous data
                    _reducedCosts.Clear();
                    var nonBasicVars = new Dictionary<string, (double min, double current, double max)>();
                    
                    for (int i = 0; i < _originalModel.Variables.Count; i++)
                    {
                        if (lpSolution.BasisIndices != null && !lpSolution.BasisIndices.Contains(i))
                        {
                            try
                            {
                                var range = _lpSensitivity.GetNonBasicVariableRange(i);
                                string varName = _originalModel.Variables[i].Name;
                                double currentValue = _solution.BestSolution[i];
                                
                                Console.WriteLine($"{varName}: [{range.lowerBound:F4}, {range.upperBound:F4}] (Current: {currentValue:F4})");
                                
                                // Store for visualization
                                nonBasicVars[varName] = (
                                    min: range.lowerBound,
                                    current: currentValue,
                                    max: range.upperBound
                                );
                                
                                // Store reduced cost for this variable
                                double reducedCost = _lpSensitivity.GetReducedCost(i);
                                _reducedCosts[varName] = reducedCost;
                            }
                            catch (Exception ex)
                            {
                                Console.WriteLine($"Error analyzing variable x{i + 1}: {ex.Message}");
                            }
                        }
                    }
                    
                    // Add visualization for non-basic variables
                    if (nonBasicVars.Count > 0)
                    {
                        Console.WriteLine("\n" + SensitivityVisualizer.CreateBarChart(
                            "NON-BASIC VARIABLES RANGES",
                            nonBasicVars,
                            width: 50
                        ));
                        
                        // Visualize non-zero reduced costs
                        var nonZeroReducedCosts = _reducedCosts
                            .Where(kvp => Math.Abs(kvp.Value) > Epsilon)
                            .ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
                            
                        if (nonZeroReducedCosts.Count > 0)
                        {
                            Console.WriteLine("\n" + SensitivityVisualizer.CreateHistogram(
                                "NON-ZERO REDUCED COSTS",
                                nonZeroReducedCosts,
                                height: 10,
                                width: 50
                            ));
                        }
                    }
                
                    // 2. Basic variable ranges
                    Console.WriteLine("\nRANGES FOR BASIC VARIABLES (AT OPTIMAL NODE):");
                    if (lpSolution.BasisIndices != null)
                    {
                        _objectiveRanges.Clear();
                        
                        foreach (int varIndex in lpSolution.BasisIndices)
                        {
                            if (varIndex >= 0 && varIndex < _originalModel.Variables.Count) // Only original variables, not slacks
                            {
                                try
                                {
                                    var range = _lpSensitivity.GetBasicVariableRange(varIndex);
                                    string varName = _originalModel.Variables[varIndex].Name;
                                    double currentValue = _solution.BestSolution[varIndex];
                                    
                                    Console.WriteLine($"{varName}: [{range.lowerBound:F4}, {range.upperBound:F4}] (Current: {currentValue:F4})");
                                    
                                    // Store for visualization
                                    _objectiveRanges[varName] = (
                                        min: range.lowerBound,
                                        current: currentValue,
                                        max: range.upperBound
                                    );
                                }
                                catch (Exception ex)
                                {
                                    Console.WriteLine($"Error analyzing basic variable x{varIndex + 1}: {ex.Message}");
                                }
                            }
                        }
                        
                        // Add visualization for objective coefficient ranges
                        if (_objectiveRanges.Count > 0)
                        {
                            Console.WriteLine("\n" + SensitivityVisualizer.CreateBarChart(
                                "OBJECTIVE COEFFICIENT RANGES",
                                _objectiveRanges,
                                width: 50
                            ));
                        }
                    }
                }

                // 3. Constraint RHS ranges
                if (_lpSensitivity != null && _solution?.CanonicalForm?.ConstraintTypes != null)
                {
                    Console.WriteLine("\nCONSTRAINT RHS RANGES (AT OPTIMAL NODE):");
                    _rhsRanges.Clear();
                    _shadowPrices.Clear();
                    
                    for (int i = 0; i < _solution.CanonicalForm.ConstraintTypes.Length; i++)
                    {
                        try
                        {
                            var range = _lpSensitivity.GetRHSRange(i);
                            double currentRHS = _solution.CanonicalForm.RightHandSide[i];
                            string constraintName = $"Constraint {i + 1} ({_solution.CanonicalForm.ConstraintTypes[i]})";
                            
                            Console.WriteLine($"{constraintName}: [{range.lowerBound:F4}, {range.upperBound:F4}] (Current: {currentRHS:F4})");
                            
                            // Store for visualization
                            _rhsRanges[constraintName] = (
                                min: range.lowerBound,
                                current: currentRHS,
                                max: range.upperBound
                            );
                            
                            // Get and store shadow price
                            try
                            {
                                double shadowPrice = _lpSensitivity.GetShadowPrice(i);
                                _shadowPrices[constraintName] = shadowPrice;
                            }
                            catch (Exception ex)
                            {
                                Console.WriteLine($"  Error getting shadow price: {ex.Message}");
                            }
                        }
                        catch (Exception ex)
                        {
                            Console.WriteLine($"Error analyzing constraint {i + 1}: {ex.Message}");
                        }
                    }
                    
                    // Add visualization for RHS ranges
                    if (_rhsRanges.Count > 0)
                    {
                        Console.WriteLine("\n" + SensitivityVisualizer.CreateBarChart(
                            "RIGHT-HAND SIDE RANGES",
                            _rhsRanges,
                            width: 50
                        ));
                    }
                    
                    // Add visualization for shadow prices
                    if (_shadowPrices.Count > 0)
                    {
                        Console.WriteLine("\n" + SensitivityVisualizer.CreateHistogram(
                            "SHADOW PRICES (DUAL VARIABLES)",
                            _shadowPrices,
                            height: 10,
                            width: 50
                        ));
                    }
                }

                // 4. Shadow prices
                if (_lpSensitivity != null && _solution?.CanonicalForm?.ConstraintTypes != null)
                {
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
                }

                // 5. Duality analysis
                if (_lpSensitivity != null)
                {
                    try
                    {
                        var dualityResult = _lpSensitivity.AnalyzeDuality();
                        if (dualityResult != null)
                        {
                            Console.WriteLine("\nDUALITY ANALYSIS:");
                            Console.WriteLine($"Primal Objective: {dualityResult.PrimalObjective:F6}");
                            Console.WriteLine($"Dual Objective: {dualityResult.DualObjective:F6}");
                            Console.WriteLine($"Duality Gap: {dualityResult.DualityGap:F6}");
                            Console.WriteLine($"Is Strong Duality: {dualityResult.HasStrongDuality}");
                        }
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine($"\nError performing duality analysis: {ex.Message}");
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error performing LP-based sensitivity analysis: {ex.Message}");
            }
        }

        private void PerformIntegerSensitivityAnalysis()
        {
            if (_originalModel?.Variables == null || _solution?.BestSolution == null)
            {
                Console.WriteLine("Cannot perform integer sensitivity analysis: Missing required data.");
                return;
            }

            try
            {
                // 1. Integer variable sensitivity
                Console.WriteLine("\nINTEGER VARIABLE SENSITIVITY:");
                for (int i = 0; i < _originalModel.Variables.Count; i++)
                {
                    if (i >= _solution.BestSolution.Length)
                    {
                        Console.WriteLine($"Warning: Variable index {i} is out of range for the solution vector.");
                        continue;
                    }

                    var variable = _originalModel.Variables[i];
                    if (variable == null)
                    {
                        Console.WriteLine($"Warning: Variable at index {i} is null.");
                        continue;
                    }

                    if (variable.Type == VariableType.Integer || variable.Type == VariableType.Binary)
                    {
                        try
                        {
                            double value = _solution.BestSolution[i];
                            double floor = Math.Floor(value);
                            double ceil = Math.Ceiling(value);
                            
                            Console.WriteLine($"x{i + 1}: Value = {value:F0} (Integer)");
                            
                            // Check if we can increase/decrease by 1 and stay integer feasible
                            bool canIncrease = true;
                            bool canDecrease = true;
                            
                            // For binary variables, we can't go below 0 or above 1
                            if (variable.Type == VariableType.Binary)
                            {
                                canIncrease = value < 1.0 - Epsilon;
                                canDecrease = value > Epsilon;
                                
                                Console.WriteLine("  Binary variable - " + 
                                    $"Can {(canIncrease ? "" : "not ")}increase to 1, " +
                                    $"Can {(canDecrease ? "" : "not ")}decrease to 0");
                            }
                            else
                            {
                                Console.WriteLine("  Integer variable - " + 
                                    $"Can increase to {value + 1}, " +
                                    $"Can decrease to {Math.Max(0, value - 1)}");
                            }
                        }
                        catch (Exception ex)
                        {
                            Console.WriteLine($"Error analyzing variable x{i + 1}: {ex.Message}");
                        }
                    }
                }

                // 2. Constraint tightness analysis
                Console.WriteLine("\nCONSTRAINT TIGHTNESS ANALYSIS:");
                var activeConstraints = new List<int>();
                var slackValues = new Dictionary<int, double>();
                
                if (_solution.CanonicalForm == null || 
                    _solution.CanonicalForm.ConstraintTypes == null || 
                    _solution.CanonicalForm.CoefficientMatrix == null ||
                    _solution.CanonicalForm.RHSVector == null)
                {
                    Console.WriteLine("  Cannot analyze constraint tightness: Missing constraint data.");
                    return;
                }
                
                for (int i = 0; i < _solution.CanonicalForm.ConstraintTypes.Length; i++)
                {
                    try
                    {
                        if (i >= _solution.CanonicalForm.CoefficientMatrix.Length || 
                            i >= _solution.CanonicalForm.RHSVector.Length)
                        {
                            Console.WriteLine($"  Warning: Constraint index {i} is out of range for coefficient matrix or RHS vector.");
                            continue;
                        }
                        
                        var coefficients = _solution.CanonicalForm.CoefficientMatrix[i];
                        if (coefficients == null)
                        {
                            Console.WriteLine($"  Warning: Coefficient matrix row {i} is null.");
                            continue;
                        }
                        
                        double lhs = 0;
                        for (int j = 0; j < Math.Min(coefficients.Length, _solution.BestSolution.Length); j++)
                        {
                            lhs += coefficients[j] * _solution.BestSolution[j];
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
                            default:
                                Console.WriteLine($"  Warning: Unknown constraint type at index {i}.");
                                continue;
                        }
                        
                        slackValues[i] = slack;
                        
                        if (slack < Epsilon)
                        {
                            activeConstraints.Add(i);
                            Console.WriteLine($"  Constraint {i + 1}: Active (Slack = {slack:E2})");
                        }
                        else
                        {
                            Console.WriteLine($"  Constraint {i + 1}: Inactive (Slack = {slack:E2})");
                        }
                    }
                    catch (Exception ex)
                    {
                        }
                }

                // 3. Branching variable impact analysis
                if (_solution.BestCandidate?.BranchingConstraint != null)
                {
                    try
                    {
                        int varIndex = -1;
                        var branchStr = _solution.BestCandidate.BranchingConstraint;
                        var match = System.Text.RegularExpressions.Regex.Match(branchStr ?? string.Empty, @"x(\d+)");
                        if (match.Success && int.TryParse(match.Groups[1].Value, out var parsedIdx))
                        {
                            varIndex = parsedIdx - 1; // convert to zero-based index
                        }
                        
                        // Validate variable index
                        if (varIndex < 0 || varIndex >= _originalModel.Variables.Count)
                        {
                            Console.WriteLine("\nBRANCHING VARIABLE ANALYSIS:");
                            Console.WriteLine($"Warning: Invalid variable index {varIndex} in branching constraint.");
                        }
                        else if (varIndex >= _solution.BestSolution.Length)
                        {
                            Console.WriteLine("\nBRANCHING VARIABLE ANALYSIS:");
                            Console.WriteLine($"Warning: Variable index {varIndex} is out of range for the solution vector.");
                        }
                        else
                        {
                            var variable = _originalModel.Variables[varIndex];
                            double value = _solution.BestSolution[varIndex];
                            double integerPart = Math.Floor(value);
                            double fractionalPart = value - integerPart;

                            Console.WriteLine("\nBRANCHING VARIABLE ANALYSIS:");
                            Console.WriteLine($"Variable x{varIndex + 1} was branched on with value: {value:F6}");
                            Console.WriteLine($"  Integer part: {integerPart}");
                            Console.WriteLine($"  Fractional part: {fractionalPart:F6}");
                            
                            if (variable.Type == VariableType.Binary)
                            {
                                Console.WriteLine("  This is a binary variable that was fixed during branching.");
                            }
                            else if (variable.Type == VariableType.Integer)
                            {
                                Console.WriteLine("  This is an integer variable that was branched on.");
                            }
                            else
                            {
                                Console.WriteLine("  Warning: This is a continuous variable that shouldn't be branched on.");
                            }
                            
                            // Analyze branching impact if we have access to the solution tree
                            if (_solution.AllNodes != null && _solution.BestCandidate.ParentId >= 0)
                            {
                                var parentNode = _solution.AllNodes
                                    .FirstOrDefault(n => n.Id == _solution.BestCandidate.ParentId);
                                    
                                if (parentNode != null)
                                {
                                    var siblings = _solution.AllNodes
                                        .Where(n => n.ParentId == parentNode.Id && n.Id != _solution.BestCandidate.Id)
                                        .ToList();
                                        
                                    if (siblings.Any())
                                    {
                                        Console.WriteLine("\nBRANCHING IMPACT ANALYSIS:");
                                        Console.WriteLine($"Found {siblings.Count} sibling nodes to compare with.");
                                        
                                        // Compare objective values with siblings
                                        var siblingObjectives = siblings
                                            .Where(s => s.Solution != null)
                                            .Select(s => s.Solution.ObjectiveValue)
                                            .ToList();
                                            
                                        if (siblingObjectives.Any())
                                        {
                                            double minObj = siblingObjectives.Min();
                                            double maxObj = siblingObjectives.Max();
                                            double avgObj = siblingObjectives.Average();
                                            
                                            Console.WriteLine($"Sibling nodes objective values: Min={minObj:F4}, " +
                                                $"Max={maxObj:F4}, Avg={avgObj:F4}");
                                                
                                            if (_solution.BestCandidate.Solution != null)
                                            {
                                                double currentObj = _solution.BestCandidate.Solution.ObjectiveValue;
                                                double bestSiblingObj = _solution.CanonicalForm.OptimizationType == OptimizationType.Maximize ? 
                                                    siblingObjectives.Max() : siblingObjectives.Min();
                                                
                                                Console.WriteLine($"Current node objective: {currentObj:F4}");
                                                Console.WriteLine($"Best sibling objective: {bestSiblingObj:F4}");
                                                
                                                // Calculate improvement over best sibling
                                                double improvement = _solution.CanonicalForm.OptimizationType == OptimizationType.Maximize ?
                                                    (currentObj - bestSiblingObj) : (bestSiblingObj - currentObj);
                                                    
                                                if (improvement > 0)
                                                {
                                                    Console.WriteLine($"Branching on x{varIndex + 1} improved the solution by {improvement:F6}");
                                                }
                                                else if (improvement < 0)
                                                {
                                                    Console.WriteLine($"Warning: Branching on x{varIndex + 1} worsened the solution by {-improvement:F6}");
                                                }
                                                else
                                                {
                                                    Console.WriteLine("No significant change in objective after branching.");
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine($"Error in branching analysis: {ex.Message}");
                    }
                }
                
                // Find siblings to compare impact of branching
                if (_solution.AllNodes != null && _solution.BestCandidate?.ParentId >= 0)
                {
                    var parentNode = _solution.AllNodes.FirstOrDefault(n => n.Id == _solution.BestCandidate?.ParentId);
                    if (parentNode != null && _solution.BestCandidate != null)
                    {
                        var siblings = _solution.AllNodes
                            .Where(n => n.ParentId == parentNode.Id && n.Id != _solution.BestCandidate?.Id)
                            .ToList();
                        
                        if (siblings.Any())
                        {
                            Console.WriteLine("\nBRANCHING IMPACT ANALYSIS:");
                            Console.WriteLine($"Found {siblings.Count} sibling nodes to compare with.");
                            
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

        // <summary>
        // Helper method to find basis variables from the final tableau with improved numerical stability
        // </summary>
        // <param name="tableau">The final simplex tableau</param>
        // <param name="numOriginalVars">Number of original decision variables</param>
        // <returns>Array of column indices that form the basis</returns>
        private int[] FindBasis(double[,] tableau, int numOriginalVars)
        {
            if (tableau == null)
                throw new ArgumentNullException(nameof(tableau));
                
            if (numOriginalVars <= 0)
                throw new ArgumentOutOfRangeException(nameof(numOriginalVars), "Number of variables must be positive");
                
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            int numConstraints = rows - 1; // Exclude objective row
            
            if (rows <= 1) // Need at least objective row + one constraint
                throw new ArgumentException("Tableau must have at least two rows", nameof(tableau));
                
            if (cols <= numOriginalVars)
                throw new ArgumentException("Tableau must have more columns than the number of original variables", nameof(tableau));
            
            var basis = new List<int>();
            var usedRows = new bool[numConstraints]; // Track which constraint rows we've used (skip objective row)
            
            // First, find all columns that are part of the basis by looking for unit vectors
            for (int j = 0; j < cols - 1; j++) // Exclude RHS column
            {
                int basicRow = -1;
                bool isBasic = true;
                
                // Check if this column is a unit vector (with some tolerance)
                for (int i = 1; i < rows; i++) // Skip objective row
                {
                    if (usedRows[i - 1]) continue; // Skip rows already assigned to basis
                    
                    double absVal = Math.Abs(tableau[i, j]);
                    
                    // Check if this could be the '1' in a unit vector
                    if (absVal.AlmostEquals(1.0, EPSILON * 10)) // Slightly more tolerance for 1.0
                    {
                        if (basicRow == -1) 
                        {
                            basicRow = i - 1; // Convert to 0-based constraint index
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
                        // Non-zero value found in a non-basic column - not a unit vector
                        isBasic = false;
                        break;
                    }
                }
                
                // If this is a basic column and we found a unique row for it
                if (isBasic && basicRow != -1)
                {
                    // Mark this row as used and add the column to the basis
                    usedRows[basicRow] = true;
                    basis.Add(j);
                    
                    // If we've found all basis variables, we can stop
                    if (basis.Count == numConstraints)
                        break;
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
        
        #region ISensitivityAnalysis Implementation
        
        public override (double lowerBound, double upperBound) GetRHSRange(int constraintIndex)
        {
            if (_lpSensitivity == null)
                throw new InvalidOperationException("LP sensitivity analysis not available");
                
            if (constraintIndex < 0 || constraintIndex >= _constraintTypes.Length)
                throw new ArgumentOutOfRangeException(nameof(constraintIndex), "Constraint index out of range");
                
            return _lpSensitivity.GetRHSRange(constraintIndex);
        }
        
        public override (double lowerBound, double upperBound) GetObjectiveCoefficientRange(int variableIndex)
        {
            if (_lpSensitivity == null)
                throw new InvalidOperationException("LP sensitivity analysis not available");
                
            if (variableIndex < 0 || variableIndex >= _variableNames.Length)
                throw new ArgumentOutOfRangeException(nameof(variableIndex), "Variable index out of range");
                
            return _lpSensitivity.GetObjectiveCoefficientRange(variableIndex);
        }
        
        public override double GetReducedCost(int variableIndex)
        {
            if (_lpSensitivity == null)
                throw new InvalidOperationException("LP sensitivity analysis not available");
                
            if (variableIndex < 0 || variableIndex >= _variableNames.Length)
                throw new ArgumentOutOfRangeException(nameof(variableIndex), "Variable index out of range");
                
            return _lpSensitivity.GetReducedCost(variableIndex);
        }
        
        public override double GetShadowPrice(int constraintIndex)
        {
            if (_lpSensitivity == null)
                throw new InvalidOperationException("LP sensitivity analysis not available");
                
            if (constraintIndex < 0 || constraintIndex >= _constraintTypes.Length)
                throw new ArgumentOutOfRangeException(nameof(constraintIndex), "Constraint index out of range");
                
            return _lpSensitivity.GetShadowPrice(constraintIndex);
        }
        
        public override Dictionary<string, double> GetSolutionValues()
        {
            var solution = new Dictionary<string, double>();
            
            if (_solution?.BestSolution == null || _originalModel?.Variables == null)
                return solution;
                
            for (int i = 0; i < Math.Min(_solution.BestSolution.Length, _originalModel.Variables.Count); i++)
            {
                solution[_originalModel.Variables[i].Name] = _solution.BestSolution[i];
            }
            
            return solution;
        }
        
        #endregion
    }
}
