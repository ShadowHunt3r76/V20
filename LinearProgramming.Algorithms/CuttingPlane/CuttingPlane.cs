using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using LinearProgramming.Parsing;
using static LinearProgramming.Algorithms.PrimalSimplex.MatrixUtils;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.OutputFormatter;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using SensitivityAnalysis = LinearProgramming.Algorithms.SensitivityAnalysis.SensitivityAnalysisImpl;

namespace LinearProgramming.Algorithms.CuttingPlane
{
    public class CuttingPlane
    {
        public const double EPSILON = 1e-10;
        private static readonly string[] InequalitySymbols = new[] { "<=", ">=", "=" };
        
        private void DisplayCanonicalForm(LinearProgramming.Parsing.CanonicalLinearProgrammingModel model)
        {
            _output.AppendLine("\nCANONICAL FORM:");
            _output.AppendLine("Objective: " + (model.OptimizationType == OptimizationType.Maximize ? "Maximize" : "Minimize") + " z = " + 
                           string.Join(" + ", model.ObjectiveCoefficients.Select((c, i) => $"{c:F2}x{i + 1}")));
            
            _output.AppendLine("\nSubject to:");
            for (int i = 0; i < model.CoefficientMatrix.Length; i++)
            {
                _output.AppendLine(FormatConstraint(model.CoefficientMatrix[i], model.RightHandSide[i], model.ConstraintTypes[i]));
            }
            
            _output.AppendLine("\nVariable Bounds:");
            for (int i = 0; i < model.ObjectiveCoefficients.Length; i++)
            {
                _output.AppendLine($"x{i + 1} >= 0");
            }
        }
        
        private string FormatConstraint(double[] constraint, double rhs, ConstraintType type)
        {
            if (constraint == null) return string.Empty;
            
            var sb = new StringBuilder();
            bool first = true;
            
            for (int i = 0; i < constraint.Length; i++)
            {
                double coef = constraint[i];
                if (Math.Abs(coef) < EPSILON) continue;
                
                if (!first && coef > 0)
                    sb.Append(" + ");
                else if (coef < 0)
                    sb.Append(" - ");
                
                if (Math.Abs(Math.Abs(coef) - 1) > EPSILON || i == constraint.Length - 1)
                    sb.Append($"{Math.Abs(coef):F2}");
                
                sb.Append($"x{i + 1}");
                first = false;
            }
            
            if (first)
                sb.Append("0");
                
            sb.Append($" {GetInequalitySymbol(type)} {rhs:F2}");
            return sb.ToString();
        }
        
        private string GetInequalitySymbol(ConstraintType type)
        {
            return type switch
            {
                ConstraintType.LessThanOrEqual => "<=",
                ConstraintType.GreaterThanOrEqual => ">=",
                ConstraintType.Equal => "=",
                _ => "?"
            };
        }
        private int _numOriginalConstraints;
        private readonly StringBuilder _output = new StringBuilder();
        private int _iterationCount = 0;
        private List<double[,]> CuttingPlaneHistory = new List<double[,]>();

        /// <summary>
        /// Checks if the model has any trivially infeasible constraints
        /// </summary>
        private string CheckModelInfeasibility(LinearProgramming.Parsing.CanonicalLinearProgrammingModel model)
        {
            var analysis = new StringBuilder();
            
            for (int i = 0; i < model.ConstraintTypes.Length; i++)
            {
                bool allNonPositive = true;
                bool allNonNegative = true;
                bool allZero = true;
                bool hasPositiveRHS = model.RHSVector[i] > EPSILON;
                bool hasNegativeRHS = model.RHSVector[i] < -EPSILON;

                foreach (var coef in model.CoefficientMatrix[i])
                {
                    if (coef > EPSILON)
                    {
                        allNonPositive = false;
                        allZero = false;
                    }
                    else if (coef < -EPSILON)
                    {
                        allNonNegative = false;
                        allZero = false;
                    }
                }

                // Check for trivially infeasible constraints
                if ((model.ConstraintTypes[i] == ConstraintType.GreaterThanOrEqual && allNonPositive && hasPositiveRHS) ||
                    (model.ConstraintTypes[i] == ConstraintType.LessThanOrEqual && allNonNegative && hasNegativeRHS) ||
                    (model.ConstraintTypes[i] == ConstraintType.Equal && allZero && Math.Abs(model.RHSVector[i]) > EPSILON))
                {
                    analysis.AppendLine($"- Constraint {i + 1} is trivially infeasible");
                    analysis.AppendLine(FormatConstraint(model.CoefficientMatrix[i].ToArray(), model.RHSVector[i], model.ConstraintTypes[i]));
                }
            }

            return analysis.Length > 0 ? analysis.ToString() : null;
        }

        /// <summary>
        /// Performs sensitivity analysis on the Cutting Plane solution
        /// </summary>
        /// <param name="lpRelaxationSolution">Solution to the LP relaxation</param>
        /// <param name="finalSolution">Final integer solution</param>
        /// <param name="model">Original model</param>
        public void PerformSensitivityAnalysis(
            LinearProgramming.Parsing.LinearProgramSolution lpRelaxationSolution,
            LinearProgramming.Parsing.LinearProgramSolution finalSolution,
            CanonicalLinearProgrammingModel model)
        {
            try
            {
                var sensitivityAnalysis = new CuttingPlaneSensitivityAnalysis(
                    lpRelaxationSolution,
                    finalSolution,
                    model,
                    CuttingPlaneHistory);
                
                sensitivityAnalysis.PerformAnalysis();
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\n⚠ Error performing sensitivity analysis: {ex.Message}");
            }
        }

        private ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel ConvertToInnerModel(LinearProgramming.Parsing.CanonicalLinearProgrammingModel model)
        {
            return new ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = model.CoefficientMatrix,
                RHSVector = model.RHSVector,
                ObjectiveCoefficients = model.ObjectiveCoefficients,
                VariableTypes = model.VariableTypes,
                ConstraintTypes = model.ConstraintTypes,
                VariableNames = model.VariableNames,
                OptimizationType = model.OptimizationType
            };
        }

        public LinearProgramming.Parsing.LinearProgramSolution CuttingPlaneSolve(LinearProgramming.Parsing.CanonicalLinearProgrammingModel model)
        {
            _output.Clear();
            _iterationCount = 0;
            CuttingPlaneHistory.Clear();
            
            // Check for trivially infeasible constraints
            string infeasibilityAnalysis = CheckModelInfeasibility(model);
            if (infeasibilityAnalysis != null)
            {
                _output.AppendLine(OutputFormatter.CreateHeader("INFEASIBILITY DETECTED"));
                _output.AppendLine("The problem is infeasible due to the following constraints:");
                _output.AppendLine(infeasibilityAnalysis);
                
                return new LinearProgramming.Parsing.LinearProgramSolution
                {
                    Status = "Infeasible",
                    SolutionVector = null,
                    ObjectiveValue = double.NaN,
                    TableHistory = new List<double[,]>()
                };
            }
            
            // Display initial canonical form
            DisplayCanonicalForm(model);
            
            // Solve LP relaxation using primal simplex
            _output.AppendLine(OutputFormatter.CreateHeader("SOLVING LP RELAXATION"));
            var simplex = new PrimalSimplexSolver();
            var innerModel = ConvertToInnerModel(model);
            var solution = simplex.Solve(innerModel);
            
            // Check for infeasibility or unboundedness in LP relaxation
            if (solution.Status == "Infeasible")
            {
                _output.AppendLine("\n✗ LP relaxation is infeasible. The integer problem is also infeasible.");
                return ConvertToParsingSolution(solution);
            }
            else if (solution.Status == "Unbounded")
            {
                _output.AppendLine("\n⚠ LP relaxation is unbounded. The integer problem may be unbounded or infeasible.");
                _output.AppendLine("  Try adding bounds on variables or additional constraints.");
                return ConvertToParsingSolution(solution);
            }
            
            double[,] tableau = solution.OptimalTable;
            
            // Display initial solution
            _output.AppendLine(OutputFormatter.CreateHeader("INITIAL LP SOLUTION"));
            _output.AppendLine($"Status: {solution.Status}");
            _output.AppendLine($"Objective Value: {solution.ObjectiveValue:F6}");
            _output.AppendLine("Solution Vector:");
            for (int i = 0; i < solution.SolutionVector.Length; i++)
            {
                _output.AppendLine(OutputFormatter.FormatKeyValue($"x{i + 1}", solution.SolutionVector[i]));
            }
            _output.AppendLine();

            // how many original constraints was started with
            _numOriginalConstraints = tableau.GetLength(0) - 1;

            int numOriginalVars = model.ObjectiveCoefficients.Length;

            // Save initial tableau (LP relaxation)
            CuttingPlaneHistory.Add((double[,])tableau.Clone());

            bool integerSolutionFound = false;
            const int maxIter = 250;
            int iter = 0;

            //Cutting plane loop
            _output.AppendLine(OutputFormatter.CreateHeader("CUTTING PLANE ITERATIONS"));
            while (!integerSolutionFound && iter < maxIter)
            {
                _iterationCount++;
                _output.AppendLine(OutputFormatter.CreateHeader($"ITERATION {_iterationCount}", 80));
                
                // Check for numerical stability
                if (double.IsInfinity(tableau[0, tableau.GetLength(1) - 1]) || 
                    double.IsNaN(tableau[0, tableau.GetLength(1) - 1]))
                {
                    _output.AppendLine("⚠ Numerical instability detected in the solution.");
                    _output.AppendLine("  This may be due to:");
                    _output.AppendLine("  1. Poorly scaled coefficients in the model");
                    _output.AppendLine("  2. Near-singular basis matrices");
                    _output.AppendLine("  3. Numerical precision issues with the current tolerance settings");
                    break;
                }
                
                // Display current tableau
                _output.AppendLine("Current Tableau:");
                _output.AppendLine(FormatTableau(tableau));
                
                // Display basis and B⁻¹
                DisplayBasisAndInverse(tableau);
                
                // Check for near-integer solution
                var currentSolution = ExtractSolution(tableau, numOriginalVars);
                
                // Create Variable array from model properties
                var variables = model.ObjectiveCoefficients
                    .Select((coeff, index) => new Variable 
                    { 
                        Name = $"x{index + 1}",
                        Coefficient = coeff
                    })
                    .ToArray();
                    
                if (IsNearlyIntegerFeasible(currentSolution, variables, out int fracVar, out double fracValue))
                {
                    _output.AppendLine($"⚠ Near-integer solution detected in x{fracVar + 1} (fractional part: {fracValue:F6})");
                    _output.AppendLine("  Consider adjusting the tolerance levels or using a different algorithm.");
                }
                
                //Find fractional row among decision variables
                int cutRow = SelectPivotRowCuttingP(tableau, numOriginalVars);
                if (cutRow == -1)
                {
                    _output.AppendLine("✓ No more fractional solutions found. Integer solution obtained.");
                    integerSolutionFound = true;
                    break;
                }
                
                _output.AppendLine($"Selected row for cut: {cutRow} (fractional variable)");

                //Apply cut
                _output.AppendLine("\nApplying Gomory's fractional cut...");
                tableau = ApplyCuttingP(tableau, cutRow);
                CuttingPlaneHistory.Add((double[,])tableau.Clone());
                
                _output.AppendLine("Tableau after adding cut:");
                _output.AppendLine(FormatTableau(tableau));

                //Restore feasibility with Dual Simplex
                _output.AppendLine("\nRestoring feasibility with Dual Simplex...");
                bool feasible = false;
                int dualIter = 0;
                int maxDualIter = 200;
                
                while (!feasible && dualIter < maxDualIter)
                {
                    _output.AppendLine($"\nDual Simplex Iteration {dualIter + 1}:");
                    
                    // Check for cycling or stalling
                    if (dualIter > maxDualIter / 2)
                    {
                        double currentObj = tableau[0, tableau.GetLength(1) - 1];
                        if (dualIter % 10 == 0)
                        {
                            _output.AppendLine($"⚠ Slow progress detected after {dualIter} iterations");
                            _output.AppendLine($"  Current objective: {currentObj:F6}");
                        }
                    }
                    
                    bool wasFeasible = feasible;
                    feasible = OptimalityCheckDualS(ref tableau);
                    
                    if (!feasible)
                    {
                        _output.AppendLine("Tableau after pivot:");
                        _output.AppendLine(FormatTableau(tableau));
                        _output.AppendLine("Price out vector: " + 
                            string.Join(", ", GetPriceOutVector(tableau)));
                        
                        // Check for cycling (same tableau seen before)
                        if (CuttingPlaneHistory.Count > 2 && 
                            IsTableauEqual(tableau, CuttingPlaneHistory[CuttingPlaneHistory.Count - 2]))
                        {
                            _output.AppendLine("⚠ Possible cycling detected - same tableau encountered twice");
                            _output.AppendLine("  Consider using a different pivot selection rule or perturbation.");
                            break;
                        }
                    }
                    
                    CuttingPlaneHistory.Add((double[,])tableau.Clone());
                    dualIter++;
                    
                    // Emergency break if we're not making progress
                    if (dualIter >= maxDualIter)
                    {
                        _output.AppendLine("⚠ Maximum dual simplex iterations reached. The problem may be ill-conditioned.");
                        break;
                    }
                }
                
                if (dualIter >= 200)
                {
                    _output.AppendLine("Warning: Dual simplex iteration limit reached.");
                }
                else
                {
                    _output.AppendLine("Feasibility restored.");
                }

                iter++;
            }

            // Display final results
            _output.AppendLine(OutputFormatter.CreateHeader("FINAL RESULTS"));
            
            if (integerSolutionFound)
            {
                _output.AppendLine("✓ Integer optimal solution found!");
                var solutionVector = ExtractSolution(tableau, numOriginalVars);
                _output.AppendLine($"Objective Value: {tableau[0, tableau.GetLength(1) - 1]:F6}");
                _output.AppendLine("Solution Vector:");
                for (int i = 0; i < solutionVector.Length; i++)
                {
                    _output.AppendLine(OutputFormatter.FormatKeyValue($"x{i + 1}", solutionVector[i]));
                }
            }
            else
            {
                _output.AppendLine("✗ Maximum iterations reached without finding an integer solution.");
            }

            // Print to console and save to file
            Console.WriteLine(_output.ToString());
            // Print final results
            PrintAnswersToTxt("cuttingplane_output.txt");

            // Perform sensitivity analysis if we have an optimal solution
            if (integerSolutionFound)
            {
                try
                {
                    _output.AppendLine("\n" + OutputFormatter.CreateHeader("PERFORMING SENSITIVITY ANALYSIS"));
                    
                    // Create Variable array from objective coefficients
                    var variables = model.ObjectiveCoefficients
                        .Select((coeff, index) => new Variable 
                        { 
                            Name = $"x{index + 1}", 
                            Coefficient = coeff,
                            LowerBound = 0,  // Assuming non-negative variables
                            UpperBound = double.PositiveInfinity
                        })
                        .ToArray();
                    
                    PerformSensitivityAnalysis(tableau, variables, model.ConstraintTypes, 
                                            numOriginalVars, _numOriginalConstraints);
                }
                catch (Exception ex)
                {
                    _output.AppendLine($"\n⚠ Error performing sensitivity analysis: {ex.Message}");
                    _output.AppendLine("  The solution is still valid, but sensitivity analysis is not available.");
                }
            }

            var finalSolution = new LinearProgramming.Parsing.LinearProgramSolution
            {
                Status = integerSolutionFound ? "Optimal Integer Solution (Cutting Plane)" : "Max iterations reached",
                SolutionVector = ExtractSolution(tableau, numOriginalVars),
                ObjectiveValue = tableau[0, tableau.GetLength(1) - 1],
                VariableValues = new Dictionary<string, double>(),
                OptimalTable = tableau,
                TableHistory = CuttingPlaneHistory
            };
            
            return finalSolution;
        }

        /// <summary>
        /// Checks if two tableaus are numerically equal within tolerance
        /// </summary>
        private void DisplayBasisAndInverse(double[,] tableau)
        {
            try
            {
                int m = tableau.GetLength(0) - 1; // Number of constraints
                int n = tableau.GetLength(1) - 1; // Number of variables including RHS
                
                _output.AppendLine("\nCurrent Basis:");
                
                // Find basic variables (columns that are part of the identity matrix)
                var basis = new List<int>();
                for (int col = 0; col < n; col++)
                {
                    int ones = 0;
                    int oneRow = -1;
                    bool isBasic = true;
                    
                    for (int row = 1; row <= m; row++)
                    {
                        if (Math.Abs(tableau[row, col] - 1) < EPSILON)
                        {
                            ones++;
                            oneRow = row;
                        }
                        else if (Math.Abs(tableau[row, col]) > EPSILON)
                        {
                            isBasic = false;
                            break;
                        }
                    }
                    
                    if (isBasic && ones == 1)
                    {
                        basis.Add(col);
                        _output.AppendLine($"  x{col + 1} = {tableau[oneRow, n]:F4}");
                    }
                }
                
                _output.AppendLine("Basis Indices: [" + string.Join(", ", basis.Select(i => $"x{i + 1}")) + "]");
            }
            catch (Exception ex)
            {
                _output.AppendLine($"\n⚠ Error displaying basis information: {ex.Message}");
            }
        }
        
        private bool IsTableauEqual(double[,] tableau1, double[,] tableau2)
        {
            if (tableau1.GetLength(0) != tableau2.GetLength(0) || 
                tableau1.GetLength(1) != tableau2.GetLength(1))
            {
                return false;
            }

            for (int i = 0; i < tableau1.GetLength(0); i++)
            {
                for (int j = 0; j < tableau1.GetLength(1); j++)
                {
                    if (Math.Abs(tableau1[i, j] - tableau2[i, j]) > EPSILON)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        
        /// <summary>
        /// Gets the price out vector (reduced costs) from the tableau
        /// </summary>
        /// <param name="tableau">The simplex tableau</param>
        /// <returns>An array of reduced costs for each non-basic variable</returns>
        private double[] GetPriceOutVector(double[,] tableau)
        {
            int numVars = tableau.GetLength(1) - 1; // Exclude RHS column
            double[] priceOutVector = new double[numVars];
            
            // The price out vector is the top row of the tableau (objective row)
            // excluding the RHS value
            for (int j = 0; j < numVars; j++)
            {
                priceOutVector[j] = tableau[0, j];
            }
            
            return priceOutVector;
        }

        /// <summary>
        /// Converts a PrimalSimplex.LinearProgramSolution to a Parsing.LinearProgramSolution
        /// </summary>
        private LinearProgramming.Parsing.LinearProgramSolution ConvertToParsingSolution(LinearProgramming.Algorithms.PrimalSimplex.LinearProgramSolution solution)
        {
            if (solution == null)
                return null;

            return new LinearProgramming.Parsing.LinearProgramSolution
            {
                Status = solution.Status,
                ObjectiveValue = solution.ObjectiveValue,
                SolutionVector = solution.SolutionVector,
                VariableValues = solution.VariableValues ?? new Dictionary<string, double>(),
                OptimalTable = solution.OptimalTable,
                TableHistory = solution.TableHistory
            };
        }
        
        /// <summary>
        /// Checks if a solution is nearly integer feasible
        /// </summary>
        private bool IsNearlyIntegerFeasible(double[] solution, Variable[] variables, out int fractionalVar, out double fractionalPart)
        {
            fractionalVar = -1;
            fractionalPart = 0.0;
            double maxFraction = 0.0;

            for (int i = 0; i < variables.Length; i++)
            {
                if (variables[i].IsInteger)
                {
                    double value = solution[i];
                    
                    // Handle very small values that might be numerical noise
                    if (Math.Abs(value) < EPSILON)
                    {
                        value = 0.0;
                    }
                    
                    double intPart = Math.Floor(value + 0.5);
                    double frac = Math.Abs(value - intPart);
                    
                    // Use relative tolerance for numbers with large magnitude
                    double tolerance = EPSILON * (1 + Math.Abs(value));
                    
                    if (frac > tolerance && frac > maxFraction)
                    {
                        maxFraction = frac;
                        fractionalVar = i;
                        fractionalPart = value - Math.Floor(value);
                    }
                }
            }

            // Consider it integer feasible if the maximum fraction is within tolerance
            return maxFraction < EPSILON * 10; 
        }

        /// <summary>
        /// Performs sensitivity analysis on the solved linear programming problem
        /// </summary>
        public void PerformSensitivityAnalysis(double[,] finalTableau, Variable[] variables, 
                                            ConstraintType[] constraintTypes, int numVariables, int numConstraints)
        {
            if (finalTableau == null || variables == null || constraintTypes == null)
                throw new ArgumentNullException("Input parameters cannot be null");

            // Determine the basis variables from the final tableau
            int[] basis = FindBasis(finalTableau, numVariables);
            
            // Create sensitivity analysis instance
            var sensitivity = new LinearProgramming.Algorithms.SensitivityAnalysis.SensitivityAnalysisImpl(
                finalTableau, basis, variables, constraintTypes, numVariables, numConstraints);
            
            _output.AppendLine("\n=== SENSITIVITY ANALYSIS ===");
            
            // 1. Display ranges for non-basic variables
            _output.AppendLine("\nRanges for Non-Basic Variables:");
            for (int i = 0; i < numVariables; i++)
            {
                if (!basis.Contains(i))
                {
                    try
                    {
                        var range = sensitivity.GetNonBasicVariableRange(i);
                        _output.AppendLine($"Variable x{i+1}: [{range.lowerBound:F4}, {range.upperBound:F4}]");
                    }
                    catch (Exception ex)
                    {
                        _output.AppendLine($"Error analyzing variable x{i+1}: {ex.Message}");
                    }
                }
            }
            
            // 2. Display ranges for basic variables
            _output.AppendLine("\nRanges for Basic Variables:");
            foreach (int varIndex in basis.Where(b => b < numVariables)) // Only original variables, not slacks
            {
                try
                {
                    var range = sensitivity.GetBasicVariableRange(varIndex);
                    _output.AppendLine($"Variable x{varIndex+1}: [{range.lowerBound:F4}, {range.upperBound:F4}]");
                }
                catch (Exception ex)
                {
                    _output.AppendLine($"Error analyzing basic variable x{varIndex+1}: {ex.Message}");
                }
            }
            
            // 3. Display ranges for constraint RHS
            _output.AppendLine("\nRanges for Constraint RHS:");
            for (int i = 0; i < numConstraints; i++)
            {
                try
                {
                    var range = sensitivity.GetRHSRange(i);
                    _output.AppendLine($"Constraint {i+1} ({constraintTypes[i]}): [{range.lowerBound:F4}, {range.upperBound:F4}]");
                }
                catch (Exception ex)
                {
                    _output.AppendLine($"Error analyzing constraint {i+1}: {ex.Message}");
                }
            }
            
            // 4. Display shadow prices
            _output.AppendLine("\nShadow Prices (Dual Values):");
            for (int i = 0; i < numConstraints; i++)
            {
                try
                {
                    double shadowPrice = sensitivity.GetShadowPrice(i);
                    _output.AppendLine($"Constraint {i+1}: {shadowPrice:F4}");
                }
                catch (Exception ex)
                {
                    _output.AppendLine($"Error getting shadow price for constraint {i+1}: {ex.Message}");
                }
            }
            
            _output.AppendLine("\n=== END OF SENSITIVITY ANALYSIS ===");
        }
        
        /// <summary>
        /// Finds the basis variables from the final tableau
        /// </summary>
        private int[] FindBasis(double[,] tableau, int numVariables)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            var basis = new List<int>();
            
            for (int j = 0; j < numVariables; j++)
            {
                int basicRow = -1;
                bool isBasic = true;

                for (int i = 1; i < rows; i++) // Skip objective row
                {
                    if (Math.Abs(tableau[i, j] - 1.0) < 1e-6)
                    {
                        if (basicRow == -1) basicRow = i - 1; // Convert to 0-based constraint index
                        else { isBasic = false; break; }
                    }
                    else if (Math.Abs(tableau[i, j]) > 1e-6)
                    {
                        isBasic = false;
                        break;
                    }
                }

                if (isBasic && basicRow != -1)
                {
                    basis.Add(j);
                }
            }
            
            return basis.ToArray();
        }

        private double[] ExtractSolution(double[,] tableau, int numOriginalVars)
        {
            if (tableau == null || tableau.Length == 0)
                return new double[numOriginalVars];
                
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            int rhsIndex = cols - 1;

            double[] solution = new double[numOriginalVars];

            for (int j = 0; j < numOriginalVars; j++)
            {
                int basicRow = -1;
                bool isBasic = true;

                for (int i = 1; i < rows; i++) // skip objective row
                {
                    if (Math.Abs(tableau[i, j] - 1.0) < 1e-6)
                    {
                        if (basicRow == -1) basicRow = i;
                        else { isBasic = false; break; }
                    }
                    else if (Math.Abs(tableau[i, j]) > 1e-6)
                    {
                        isBasic = false;
                        break;
                    }
                }

                if (isBasic && basicRow != -1)
                {
                    solution[j] = tableau[basicRow, rhsIndex];
                    
                    // Check for numerical instability in basic variables
                    if (double.IsInfinity(solution[j]) || double.IsNaN(solution[j]))
                    {
                        _output.AppendLine($"⚠ Warning: Basic variable x{j+1} has invalid value: {solution[j]}");
                    }
                }
                else
                {
                    solution[j] = 0.0;
                }
            }

            return solution;
        }

        //cutting plane pivot row selecection
        private int SelectPivotRowCuttingP(double[,] tableau, int numOriginalVariables)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            int rhsIndex = cols - 1;

            const double tol = NumericalStabilityUtils.Epsilon; // tolerance to avoid endless tiny fractions
            int bestRow = -1;
            double bestFraction = 0.0;

            // Only scan ORIGINAL constraints, not the cuts added later
            int lastOriginalRow = Math.Min(_numOriginalConstraints, rows - 1);
            for (int i = 1; i <= lastOriginalRow; i++)
            {
                double rhs = tableau[i, rhsIndex];
                double frac = rhs - Math.Floor(rhs);
                if (frac < tol || frac > 1 - tol) continue; // treat near-integers as integers

                // find the basic column for this row (unit column test)
                int basicCol = -1;
                for (int j = 0; j < rhsIndex; j++)
                {
                    if (Math.Abs(tableau[i, j] - 1.0) < 1e-6)
                    {
                        bool isUnit = true;

                        // check other rows for 0
                        for (int k = 0; k < rows; k++)
                        {
                            if (k == i) continue;
                            if (Math.Abs(tableau[k, j]) > 1e-6) { isUnit = false; break; }
                        }
                        if (isUnit) { basicCol = j; break; }
                    }
                }

                // only cut if the basic var is one of the ORIGINAL x-variables
                if (basicCol != -1 && basicCol < numOriginalVariables)
                {
                    if (frac > bestFraction)
                    {
                        bestFraction = frac;
                        bestRow = i;
                    }
                }
            }

            return bestRow; // -1 => no fractional original x found
        }


        //cutting plane cut method

        private double[,] ApplyCuttingP(double[,] tableau, int cutRow)
        {

            // Helper: fractional part function
            double FractionalPart(double x)
            {
                double f = x - Math.Floor(x); // always between 0 and 1 
                return (Math.Abs(f) < 1e-6) ? 0.0 : f; // treat tiny value as 0 
            }

            int oldRows = tableau.GetLength(0);
            int oldCols = tableau.GetLength(1);
            int oldRhsIndex = oldCols - 1;

            // New matrix has one extra row (the cut) and one extra column (new slack). 
            // Insert new slack column before RHS, so RHS stays the last column. 
            int newRows = oldRows + 1;
            int newCols = oldCols + 1;        // one extra column
            int newRhsIndex = newCols - 1;    // RHS stays last

            double[,] newT = new double[newRows, newCols];

            // Copy old coefficients, shifting the old RHS one column to the RIGHT.
            // Columns [0 .. oldRhsIndex-1] stay the same indices.
            // The NEW slack column is at index = oldRhsIndex.
            for (int i = 0; i < oldRows; i++)
            {
                // copy all non-RHS columns as-is
                for (int j = 0; j < oldRhsIndex; j++)
                    newT[i, j] = tableau[i, j];

                // set the NEW slack column to 0 for existing rows
                newT[i, oldRhsIndex] = 0.0;

                // move the RHS from oldRhsIndex -> newRhsIndex
                newT[i, newRhsIndex] = tableau[i, oldRhsIndex];
            }

            // Build the new cut row (index = oldRows)
            // coefficients (negative fractional parts) for all decision/slack columns (excluding RHS)
            for (int j = 0; j < oldRhsIndex; j++)
                newT[oldRows, j] = -FractionalPart(tableau[cutRow, j]);

            // set the NEW slack column coefficient to +1
            newT[oldRows, oldRhsIndex] = 1.0;

            // RHS = - fractional part of the original RHS
            newT[oldRows, newRhsIndex] = -FractionalPart(tableau[cutRow, oldRhsIndex]);

            return newT;
        }



        //Dual Simplex programming

        // pivot row selection DualSimplex

        private int SelectPivotRowDualS(double[,] tableau)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            int rhsIndex = cols - 1;

            int pivotRow = -1;
            double mostNegative = 0.0; // we want < 0

            // Check constraint rows (skip row 0 = objective)
            for (int i = 1; i < rows; i++)
            {
                double rhs = tableau[i, rhsIndex];
                if (rhs < mostNegative - 1e-6) 
                {
                    mostNegative = rhs;
                    pivotRow = i;
                }
            }

            return pivotRow; // -1 means "no negative RHS found"
        }

        //pivot column selection DualSimplex

        private int SelectPivotColumnDualS(double[,] tableau, int pivotRow)
        {
            int cols = tableau.GetLength(1);
            int rhsIndex = cols - 1;

            int pivotCol = -1; //indicating pivotCol not found as starting point
            double minRatio = double.MaxValue; //starting point for "smallest ratio" is largest possible double value

            for (int j = 0; j < rhsIndex; j++)
            {
                double coeff = tableau[pivotRow, j];

                if (coeff < -1e-6) // only consider negative pivot row entries
                {
                    double ratio = tableau[0, j] / coeff;

                    if (ratio >= 0 && ratio < minRatio) // only accept non-negative ratios
                    {
                        minRatio = ratio;
                        pivotCol = j;
                    }
                }
            }

            return pivotCol; // -1 means no valid column found
        }

        // optimality check Dual Simplex

        private bool OptimalityCheckDualS(ref double[,] tableau)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            int rhsIndex = cols - 1;

            // Step 1: check feasibility (all RHS >= 0)
            for (int i = 1; i < rows; i++) // skip row 0 (objective)
            {
                if (tableau[i, rhsIndex] < -1e-6) // negative RHS found
                {
                    Console.WriteLine("Dual Simplex: Solution infeasible, selecting pivot...");

                    int pivotRow = SelectPivotRowDualS(tableau);
                    if (pivotRow == -1)
                    {
                        Console.WriteLine("No valid pivot row found. Problem infeasible.");
                        return false;
                    }

                    int pivotCol = SelectPivotColumnDualS(tableau, pivotRow);
                    if (pivotCol == -1)
                    {
                        Console.WriteLine("No valid pivot column found. Problem infeasible.");
                        return false;
                    }

                    Console.WriteLine($"Pivoting on row {pivotRow}, column {pivotCol}...");
                    tableau = Pivot(tableau, pivotCol, pivotRow);

                    // return false = "not optimal yet, keep looping"
                    return false;
                }
            }

            //Here all RHS are ≥ 0
            Console.WriteLine("Dual Simplex: Feasible solution found.");
           
            return true; // solution feasible
        }


        // Local copy of Pivot so CuttingPlane is self-contained (else we can make Pivot in PrimalSimplexSolver public and call here)
        private double[,] Pivot(double[,] tableau, int entering, int leaving)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            double[,] newTable = new double[rows, cols];
            double pivot = tableau[leaving, entering];

            
            for (int j = 0; j < cols; j++)
                newTable[leaving, j] = tableau[leaving, j] / pivot;

            // Eliminate column entries
            for (int i = 0; i < rows; i++)
            {
                if (i == leaving) continue;
                double factor = tableau[i, entering];
                for (int j = 0; j < cols; j++)
                    newTable[i, j] = tableau[i, j] - factor * newTable[leaving, j];
            }

            return newTable;
        }


        private void PrintAnswersToTxt(string filePath = "cuttingplane_output.txt")
        {
            using (StreamWriter writer = new StreamWriter(filePath))
            {
                writer.WriteLine("=== Cutting Plane Method Output ===\n");

                int iter = 1;
                foreach (var tableau in CuttingPlaneHistory)
                {
                    writer.WriteLine($"--- Cutting Plane Iteration {iter} ---");
                    PrintMatrix(writer, tableau);
                    writer.WriteLine();
                    iter++;
                }

                if (CuttingPlaneHistory.Count > 0)
                {
                    writer.WriteLine("=== Final Integer Optimal Tableau ===");
                    var finalTableau = CuttingPlaneHistory[CuttingPlaneHistory.Count - 1];
                    PrintMatrix(writer, finalTableau);
                }
            }

            Console.WriteLine($"Cutting plane output written to {filePath}");
        }

        // Helper to print any tableau with values rounded to three decimals
        private void PrintMatrix(StreamWriter writer, double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    writer.Write($"{Math.Round(matrix[i, j], 3),8} "); // 3 decimals, padded
                }
                writer.WriteLine();
            }
        }

    }
}
