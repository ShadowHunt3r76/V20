using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.PrimalSimplex;

// Use the outer CanonicalLinearProgrammingModel from LinearProgramming.Parsing
using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;
using LinearProgramSolution = LinearProgramming.Algorithms.PrimalSimplex.LinearProgramSolution;
using TableauIteration = LinearProgramming.Algorithms.PrimalSimplex.TableauIteration;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    public class PrimalSimplexSolver
    {
        private CanonicalLinearProgrammingModel model;
        private const double Epsilon = MatrixUtils.Epsilon;
        private string[] _variableNames;
        private List<int> _basisIndices = new List<int>();
        private int _slackCount;
        private int _artificialCount;
        private int _originalVarCount;

        public LinearProgramSolution Solve(CanonicalLinearProgrammingModel model)
        {
            // Store the model reference for use in other methods
            this.model = model ?? throw new ArgumentNullException(nameof(model), "Model cannot be null");
                
            // Input validation
            if (model.CoefficientMatrix == null || model.CoefficientMatrix.Length == 0)
                throw new ArgumentException("Coefficient matrix cannot be null or empty", nameof(model.CoefficientMatrix));
                
            if (model.ObjectiveCoefficients == null || model.ObjectiveCoefficients.Length == 0)
                throw new ArgumentException("Objective coefficients cannot be null or empty", nameof(model.ObjectiveCoefficients));
                
            if (model.RightHandSide == null || model.RightHandSide.Length == 0)
                throw new ArgumentException("Right-hand side values cannot be null or empty", nameof(model.RightHandSide));
                
            if (model.ConstraintTypes == null || model.ConstraintTypes.Length == 0)
                throw new ArgumentException("Constraint types cannot be null or empty", nameof(model.ConstraintTypes));
                
            if (model.CoefficientMatrix.Length != model.ConstraintTypes.Length || 
                model.CoefficientMatrix.Length != model.RightHandSide.Length)
                throw new ArgumentException("Mismatched array dimensions in the model");
                
            // Check for NaN or Infinity in coefficients
            foreach (var row in model.CoefficientMatrix)
            {
                if (row == null)
                    throw new ArgumentException("Coefficient matrix row cannot be null");
                    
                foreach (var coef in row)
                {
                    if (double.IsNaN(coef) || double.IsInfinity(coef))
                        throw new ArgumentException("Coefficients must be finite numbers");
                }
            }
            
            // Initialize variable names if not provided
            if (model.VariableNames == null || model.VariableNames.Length == 0)
            {
                var varNames = new List<string>();
                for (int i = 0; i < model.ObjectiveCoefficients.Length; i++)
                {
                    varNames.Add($"x{i + 1}");
                }
                model.VariableNames = varNames.ToArray();
            }
            
            // Initialize solution and display canonical form
            InitializeSolution(model, out var solution, out var tableau, out var basis);
            solution.VariableNames = _variableNames;
            
            // Display problem information
            Console.WriteLine("\n=== Solving Linear Program ===");
            DisplayCanonicalForm(model);
            Console.WriteLine("\n=== Starting Primal Simplex Method ===");
            
            int maxIter = 1000;
            int iter = 0;
            bool isDegenerate = false;

            while (!IsOptimal(tableau) && iter < maxIter)
            {
                int entering = SelectEnteringVariable(tableau);
                if (entering == -1)
                {
                    solution.Status = "Optimal";
                    break;
                }

                int leaving = SelectLeavingVariable(tableau, entering);
                if (leaving == -1)
                {
                    solution.Status = "Unbounded";
                    solution.OptimalTable = tableau;
                    return solution;
                }

                // Check for degeneracy (tie in minimum ratio test)
                if (IsDegeneratePivot(tableau, entering, leaving))
                {
                    isDegenerate = true;
                }

                // Record and display iteration before pivot
                var iteration = CreateIteration(iter + 1, tableau, basis, entering, leaving, "Pivoting");
                solution.Iterations.Add(iteration);
                
                // Display iteration information
                Console.WriteLine($"\n=== Iteration {iteration.Iteration} ===");
                Console.WriteLine($"Entering: {_variableNames[entering]} (x{entering + 1})");
                Console.WriteLine($"Leaving: {basis[leaving - 1]} (row {leaving})");
                Console.WriteLine("\nCurrent Tableau:");
                Console.WriteLine(TableauToString(tableau));

                // Perform pivot
                basis[leaving - 1] = _variableNames[entering];
                Pivot(tableau, entering, leaving);
                
                iter++;
            }

            FinalizeSolution(solution, tableau, basis, isDegenerate);
            
            // Display final solution
            Console.WriteLine("\n=== Solution Complete ===");
            solution.DisplaySolution();
            
            return solution;
        }

        private void InitializeSolution(CanonicalLinearProgrammingModel model, 
            out LinearProgramSolution solution, 
            out double[,] tableau, 
            out string[] basis)
        {
            int m = model.CoefficientMatrix.Length;
            _originalVarCount = model.ObjectiveCoefficients.Length;
            _slackCount = 0;
            _artificialCount = 0;

            // Count slack and artificial variables based on constraint types
            foreach (var constraintType in model.ConstraintTypes)
            {
                if (constraintType == ConstraintType.LessThanOrEqual || 
                    constraintType == ConstraintType.GreaterThanOrEqual)
                {
                    _slackCount++;
                }
                
                if (constraintType == ConstraintType.GreaterThanOrEqual || 
                    constraintType == ConstraintType.Equal)
                {
                    _artificialCount++;
                }
            }

            // Use provided variable names if available, otherwise generate them
            if (model.VariableNames != null && model.VariableNames.Length >= _originalVarCount + _slackCount + _artificialCount)
            {
                _variableNames = model.VariableNames;
            }
            else
            {
                // Initialize variable names if not provided or incomplete
                InitializeVariableNames(m);
            }
            
            // Build initial tableau
            tableau = BuildInitialTableau(model);
            
            // Initialize basis variables
            basis = InitializeBasisVariables(m);
            
            // Initialize solution
            solution = new LinearProgramSolution
            {
                CanonicalMatrix = MatrixUtils.ConvertToMatrix(model.CoefficientMatrix),
                VariableNames = _variableNames,
                BasisVariables = basis,
                InitialTable = (double[,])tableau.Clone(),
                Status = "In Progress"
            };
        }

    private void InitializeVariableNames(int numConstraints)
    {
        var varList = new List<string>();
        
        // Original variables (x1, x2, ..., xn)
        for (int i = 0; i < _originalVarCount; i++)
        {
            varList.Add($"x{i + 1}");
        }
        
        int slackIndex = 1;
        int artificialIndex = 1;
        
        // Add slack/surplus and artificial variables based on constraint types
        for (int i = 0; i < numConstraints; i++)
        {
            var constraintType = model.ConstraintTypes[i];
            
            // Add slack/surplus variables
            if (constraintType != ConstraintType.Equal)
            {
                varList.Add($"s{slackIndex++}");
            }
            
            // Add artificial variables for >= and = constraints
            if (constraintType == ConstraintType.GreaterThanOrEqual || 
                constraintType == ConstraintType.Equal)
            {
                varList.Add($"a{artificialIndex++}");
            }
        }
        
        _variableNames = varList.ToArray();
    }

    private string[] InitializeBasisVariables(int numConstraints)
    {
        var basis = new string[numConstraints];
        int slackIndex = 1;
        int artificialIndex = 1;
        
        for (int i = 0; i < numConstraints; i++)
        {
            var constraintType = model.ConstraintTypes[i];
            
            // For <= constraints, use slack variable in basis
            if (constraintType == ConstraintType.LessThanOrEqual)
            {
                basis[i] = $"s{slackIndex++}";
            }
            // For >= and = constraints, use artificial variable in basis
            else if (constraintType == ConstraintType.GreaterThanOrEqual || 
                    constraintType == ConstraintType.Equal)
            {
                // If there's a slack variable (for >=), it's not in basis
                if (constraintType == ConstraintType.GreaterThanOrEqual)
                {
                    slackIndex++; // Skip the slack variable for this constraint
                }
                basis[i] = $"a{artificialIndex++}";
            }
        }
        
        return basis;
    }

    /// <summary>
    /// Creates a new iteration object for the simplex method
    /// </summary>
    private TableauIteration CreateIteration(int iter, double[,] tableau, string[] basis, 
        int entering, int leaving, string status)
    {
        int numVars = _originalVarCount + _slackCount + _artificialCount;
        int numConstraints = tableau.GetLength(0) - 1;
        
        // Extract the current solution
        ExtractSolution(tableau, basis, out double[] solutionVector, out double objectiveValue);
        
        // Get the current basis inverse if available
        double[,] basisInverse = null;
        if (numConstraints > 0 && numVars > 0)
        {
            basisInverse = new double[numConstraints, numConstraints];
            for (int i = 0; i < numConstraints; i++)
            {
                for (int j = 0; j < numConstraints; j++)
                {
                    basisInverse[i, j] = tableau[i + 1, _originalVarCount + j];
                }
            }
        }
        
        // Get reduced costs
        double[] reducedCosts = new double[numVars];
        for (int j = 0; j < numVars; j++)
        {
            reducedCosts[j] = tableau[0, j];
        }
        
        // Get simplex multipliers (shadow prices)
        double[] simplexMultipliers = new double[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            simplexMultipliers[i] = -tableau[0, _originalVarCount + i];
        }
        
        // Get basis indices
        int[] basisIndices = new int[basis.Length];
        for (int i = 0; i < basis.Length; i++)
        {
            // Find the index of the basis variable in the variable names array
            int index = Array.IndexOf(_variableNames, basis[i]);
            basisIndices[i] = index >= 0 ? index : -1;
        }
        
        // Get non-basis indices
        var nonBasisVars = _variableNames.Except(basis).ToArray();
        int[] nonBasisIndices = new int[nonBasisVars.Length];
        for (int i = 0; i < nonBasisVars.Length; i++)
        {
            int index = Array.IndexOf(_variableNames, nonBasisVars[i]);
            nonBasisIndices[i] = index >= 0 ? index : -1;
        }
        
        // Create and return the iteration object
        return new TableauIteration
        {
            Iteration = iter,
            BasisVariables = new List<string>(basis),
            NonBasisVariables = nonBasisVars,
            NonBasicVariableIndices = nonBasisIndices,
            BasicVariableValues = GetBasicVariableValues(tableau, basis),
            PivotElement = leaving >= 0 && entering >= 0 ? tableau[leaving + 1, entering] : 0.0
        };
    }

    /// <summary>
    /// Extracts the solution vector and objective value from the tableau
    /// </summary>
    private void ExtractSolution(double[,] tableau, string[] basis, out double[] solutionVector, out double objectiveValue)
    {
        int numVars = _originalVarCount + _slackCount + _artificialCount;
        int numConstraints = tableau.GetLength(0) - 1; // Exclude objective row
        solutionVector = new double[numVars];
        
        // Initialize all variables to 0
        Array.Fill(solutionVector, 0);
        
        // Set basic variables
        for (int i = 0; i < numConstraints && i < basis.Length; i++)
        {
            string varName = basis[i];
            if (!string.IsNullOrEmpty(varName) && varName.Length > 0)
            {
                // Find the index of this variable in the solution vector
                int varIndex = Array.IndexOf(_variableNames, varName);
                if (varIndex >= 0 && varIndex < numVars)
                {
                    // The RHS value for this basic variable is in the last column of its row
                    solutionVector[varIndex] = tableau[i + 1, numVars];
                }
            }
        }
        
        // The objective value is in the bottom-right corner of the tableau (negated)
        objectiveValue = -tableau[0, numVars];
    }
    
    /// <summary>
    /// Gets the values of the basic variables from the tableau
    /// </summary>
    private double[] GetBasicVariableValues(double[,] tableau, string[] basis)
    {
        int numVars = _originalVarCount + _slackCount + _artificialCount;
        int numConstraints = tableau.GetLength(0) - 1;
        double[] values = new double[numConstraints];
        
        // The RHS values are in the last column of the tableau (excluding objective row)
        for (int i = 0; i < numConstraints; i++)
        {
            values[i] = tableau[i + 1, numVars];
        }
        
        return values;
    }
    
    private void FinalizeSolution(LinearProgramSolution solution, double[,] tableau, string[] basis, bool isDegenerate)
    {
        // Set the final tableau and status
        solution.OptimalTable = (double[,])tableau.Clone();
        solution.Status = IsOptimal(tableau) ? "Optimal" : "Infeasible";
        solution.IsDegenerate = isDegenerate;
        
        // Extract the solution from the final tableau
        ExtractSolution(tableau, basis, out double[] solutionVector, out double objectiveValue);
        solution.SolutionVector = solutionVector;
        solution.ObjectiveValue = objectiveValue;
        
        // Set variable values in the solution
        var variableValues = new Dictionary<string, double>();
        for (int i = 0; i < _variableNames.Length; i++)
        {
            if (i < solutionVector.Length)
            {
                variableValues[_variableNames[i]] = solutionVector[i];
            }
        }
        solution.VariableValues = variableValues;
        
        // Add final iteration to history
        solution.Iterations.Add(CreateIteration(
            solution.Iterations.Count + 1, 
            tableau, 
            basis, 
            -1, 
            -1, 
            solution.Status));
            
        // Set the basis variables in the solution
        solution.BasisVariables = basis.ToArray();
        
        // Set the final objective value
        solution.ObjectiveValue = objectiveValue;
        
        // Set the final status based on the solution
        if (solution.Status == "Optimal")
        {
            // Check for degeneracy
            if (isDegenerate)
            {
                solution.Status = "Degenerate Optimal";
            }
            
            // Check for alternative optima
            bool hasAlternativeOptima = false;
            int numVars = _originalVarCount + _slackCount + _artificialCount;
            for (int j = 0; j < numVars; j++)
            {
                if (Math.Abs(tableau[0, j]) < Epsilon && solutionVector[j] > Epsilon)
                {
                    hasAlternativeOptima = true;
                    break;
                }
            }
            
            if (hasAlternativeOptima)
            {
                solution.Status += " (Multiple Optimal Solutions)";
            }
        }
    }

    /// <summary>
    /// Displays the linear programming problem in canonical form
    /// </summary>
    private void DisplayCanonicalForm(CanonicalLinearProgrammingModel model)
    {
        Console.WriteLine("\n=== Canonical Form ===");
        
        // Display objective function
        Console.Write("Maximize: ");
        bool firstTerm = true;
        for (int j = 0; j < model.ObjectiveCoefficients.Length; j++)
        {
            double coef = model.ObjectiveCoefficients[j];
            if (Math.Abs(coef) < Epsilon) continue;
            
            string varName = j < model.VariableNames.Length ? model.VariableNames[j] : $"x{j + 1}";
            
            if (firstTerm)
            {
                firstTerm = false;
                Console.Write(coef < 0 ? "-" : "");
            }
            else
            {
                Console.Write(coef < 0 ? " - " : " + ");
            }
            
            double absCoef = Math.Abs(coef);
            if (Math.Abs(absCoef - 1.0) > Epsilon)
            {
                Console.Write($"{absCoef} ");
            }
            Console.Write(varName);
        }
        Console.WriteLine("\n");
        
        // Display constraints
        Console.WriteLine("Subject to:");
        for (int i = 0; i < model.CoefficientMatrix.Length; i++)
        {
            firstTerm = true;
            for (int j = 0; j < model.CoefficientMatrix[i].Length; j++)
            {
                double coef = model.CoefficientMatrix[i][j];
                if (Math.Abs(coef) < Epsilon) continue;
                
                string varName = j < model.VariableNames.Length ? model.VariableNames[j] : $"x{j + 1}";
                
                if (firstTerm)
                {
                    firstTerm = false;
                    Console.Write(coef < 0 ? "-" : "  ");
                }
                else
                {
                    Console.Write(coef < 0 ? " - " : " + ");
                }
                
                double absCoef = Math.Abs(coef);
                if (Math.Abs(absCoef - 1.0) > Epsilon)
                {
                    Console.Write($"{absCoef} ");
                }
                Console.Write(varName);
            }
            
            // Display constraint type and RHS
            string constraintType = "";
            switch (model.ConstraintTypes[i])
            {
                case ConstraintType.LessThanOrEqual:
                    constraintType = "<=";
                    break;
                case ConstraintType.GreaterThanOrEqual:
                    constraintType = ">=";
                    break;
                case ConstraintType.Equal:
                    constraintType = "=";
                    break;
            }
            
            Console.WriteLine($" {constraintType} {model.RHSVector[i]}");
        }
        
        // Display variable types
        Console.WriteLine("\nVariable types:");
        for (int j = 0; j < model.VariableTypes.Length; j++)
        {
            string varName = j < model.VariableNames.Length ? model.VariableNames[j] : $"x{j + 1}";
            Console.WriteLine($"{varName}: {model.VariableTypes[j]}");
        }
        
        Console.WriteLine("\n" + new string('=', 80) + "\n");
    }

    private double[,] BuildInitialTableau(CanonicalLinearProgrammingModel model)
    {
        int numConstraints = model.CoefficientMatrix.Length;
        int numOriginalVars = model.ObjectiveCoefficients.Length;
        int totalVars = _originalVarCount + _slackCount + _artificialCount;
        
        // Create tableau with extra columns for RHS
        var tableau = new double[numConstraints + 1, totalVars + 1];
        
        // 1. Set up objective row (negated for maximization)
        for (int j = 0; j < numOriginalVars; j++)
        {
            tableau[0, j] = -model.ObjectiveCoefficients[j];
        }
        
        // Initialize slack and artificial variable positions
        int slackIndex = 0;
        int artificialIndex = 0;
        
        // 2. Set up constraint rows
        for (int i = 0; i < numConstraints; i++)
        {
            // Original coefficients
            for (int j = 0; j < numOriginalVars; j++)
            {
                tableau[i + 1, j] = model.CoefficientMatrix[i][j];
            }
            
            var constraintType = model.ConstraintTypes[i];
            
            // Handle slack/surplus and artificial variables based on constraint type
            if (constraintType == ConstraintType.LessThanOrEqual)
            {
                // Add slack variable
                int slackPos = _originalVarCount + slackIndex;
                tableau[i + 1, slackPos] = 1;
                slackIndex++;
            }
            else if (constraintType == ConstraintType.GreaterThanOrEqual)
            {
                // Add surplus variable (negative slack)
                int slackPos = _originalVarCount + slackIndex;
                tableau[i + 1, slackPos] = -1;
                slackIndex++;
                
                // Add artificial variable for >= constraint
                int artPos = _originalVarCount + _slackCount + artificialIndex;
                tableau[i + 1, artPos] = 1;
                artificialIndex++;
                
                // Use Big M method for artificial variables in the objective
                tableau[0, artPos] = double.MaxValue / 2;
            }
            else if (constraintType == ConstraintType.Equal)
            {
                // Add artificial variable for = constraint
                int artPos = _originalVarCount + _slackCount + artificialIndex;
                tableau[i + 1, artPos] = 1;
                artificialIndex++;
                
                // Use Big M method for artificial variables in the objective
                tableau[0, artPos] = double.MaxValue / 2;
            }
            
            // Set RHS value for the constraint
            tableau[i + 1, totalVars] = model.RHSVector[i];
        }
        
        // 3. Make sure the objective coefficients of basic variables are zero
        for (int i = 0; i < numConstraints; i++)
        {
            // Find the basic variable (with coefficient 1) in this row
            for (int j = 0; j < totalVars; j++)
            {
                if (Math.Abs(tableau[i + 1, j] - 1.0) < Epsilon)
                {
                    // Found a basic variable, make its objective coefficient zero
                    double factor = tableau[0, j];
                    if (Math.Abs(factor) > Epsilon)
                    {
                        // Subtract factor * constraint row from objective row
                        for (int k = 0; k <= totalVars; k++)
                        {
                            tableau[0, k] -= factor * tableau[i + 1, k];
                        }
                    }
                    break;
                }
            }
        }
        
        return tableau;
    }

    private bool IsOptimal(double[,] tableau)
    {
        int cols = tableau.GetLength(1);
        for (int j = 0; j < cols - 1; j++)
        {
            if (tableau[0, j] < -Epsilon) // Using epsilon for floating-point comparison
                return false;
        }
        return true;
    }
    
    private int SelectEnteringVariable(double[,] tableau)
    {
        if (tableau == null)
            throw new ArgumentNullException(nameof(tableau), "Tableau cannot be null");
            
        int colCount = tableau.GetLength(1);
        int entering = -1;
        double maxReducedCost = -Epsilon;
        
        // Check for numerical stability
        for (int j = 0; j < colCount - 1; j++) // Exclude RHS column
        {
            double reducedCost = tableau[0, j];
            if (double.IsNaN(reducedCost) || double.IsInfinity(reducedCost))
                throw new InvalidOperationException("Numerical instability detected in reduced costs");
        }
        
        for (int j = 0; j < colCount - 1; j++)
        {
            double reducedCost = tableau[0, j];
            if (reducedCost < maxReducedCost)
            {
                maxReducedCost = reducedCost;
                entering = j;
            }
        }
        
        return entering;
    }

    private int SelectLeavingVariable(double[,] tableau, int entering)
    {
        if (tableau == null)
            throw new ArgumentNullException(nameof(tableau), "Tableau cannot be null");
            
        int rowCount = tableau.GetLength(0);
        int colCount = tableau.GetLength(1);
        
        if (entering < 0 || entering >= colCount - 1) // -1 to exclude RHS column
            throw new ArgumentOutOfRangeException(nameof(entering), "Entering variable index is out of range");
            
        int leaving = -1;
        double minRatio = double.MaxValue;
        bool foundPositive = false;
        
        // Check for numerical stability in the entering column
        bool allNonPositive = true;
        for (int i = 1; i < rowCount; i++)
        {
            double val = tableau[i, entering];
            if (double.IsNaN(val) || double.IsInfinity(val))
                throw new InvalidOperationException("Numerical instability detected in pivot column");
                
            if (val > Epsilon)
                allNonPositive = false;
        }
        
        if (allNonPositive)
            return -1; // Unbounded problem
        
        for (int i = 1; i < rowCount; i++)
        {
            double val = tableau[i, entering];
            if (val > Epsilon)
            {
                foundPositive = true;
                double ratio = tableau[i, colCount - 1] / val;
                if (ratio < minRatio)
                {
                    minRatio = ratio;
                    leaving = i;
                }
            }
        }
        
        if (!foundPositive)
            throw new InvalidOperationException("No positive values found in pivot column");
        
        return leaving;
    }

    private bool IsDegeneratePivot(double[,] tableau, int entering, int leaving)
    {
        if (leaving <= 0) return false;
        
        int rows = tableau.GetLength(0);
        int cols = tableau.GetLength(1);
        
        // Check if the RHS of the leaving row is close to zero (indicating degeneracy)
        if (Math.Abs(tableau[leaving, cols - 1]) < Epsilon)
            return true;
            
        // Check for alternative pivots with same ratio
        double minRatio = double.MaxValue;
        for (int i = 1; i < rows; i++)
        {
            if (tableau[i, entering] > Epsilon)
            {
                double ratio = tableau[i, cols - 1] / tableau[i, entering];
                if (ratio >= 0 && ratio < minRatio)
                    minRatio = ratio;
            }
        }
        
        int count = 0;
        for (int i = 1; i < rows; i++)
        {
            if (tableau[i, entering] > Epsilon)
            {
                double ratio = tableau[i, cols - 1] / tableau[i, entering];
                if (Math.Abs(ratio - minRatio) < Epsilon)
                    count++;
            }
        }
    
    return count > 1; // Return whether the pivot is degenerate (multiple minimum ratios)
}
        
/// <summary>
/// Generates a string representation of the current tableau
/// </summary>
    private string TableauToString(double[,] tableau)
    {
        if (tableau == null) return "Tableau is null";
        
        int rows = tableau.GetLength(0);
        int cols = tableau.GetLength(1);
    
        var sb = new StringBuilder();
    
        // Print header
        sb.Append("Basis\t");
        for (int j = 0; j < cols - 1; j++)
        {
            string varName = j < _variableNames.Length ? _variableNames[j] : $"x{j + 1}";
            sb.Append($"{varName,10}");
        }
        sb.AppendLine("|" + "RHS".PadLeft(10));
        
        // Print rows
        for (int i = 0; i < rows; i++)
        {
            // Print basis variable if available
            if (i > 0 && _basisIndices != null && i - 1 < _basisIndices.Count)
            {
                int varIndex = _basisIndices[i - 1];
                string basisVar = GetVariableName(varIndex);
                sb.Append($"{basisVar,-5}\t");
            }
            else
            {
                sb.Append("     \t");
            }
            
            // Print row values
            for (int j = 0; j < cols; j++)
            {
                if (j == cols - 1)
                    sb.Append("|");
                    
                // Format the number with fixed width and decimal places
                string valStr = Math.Abs(tableau[i, j]) < Epsilon ? "0.000" : $"{tableau[i, j],10:0.000}";
                sb.Append(valStr);
            }
            sb.AppendLine();
            
            // Print separator after objective row
            if (i == 0)
        {
            sb.Append(new string('-', 10 * (cols + 1) + 5));
            sb.AppendLine();
        }
    }
    
    return sb.ToString();
    }
    
    /// <summary>
    /// Gets the name of a variable by its index
    /// </summary>
    private string GetVariableName(int index)
    {
        if (index < _originalVarCount)
        {
            return $"x{index + 1}";
        }
        else if (index < _originalVarCount + _slackCount)
        {
            return $"s{index - _originalVarCount + 1}";
        }
        else
        {
            return $"a{index - _originalVarCount - _slackCount + 1}";
        }
    }
    
    /// <summary>
    /// Performs the pivot operation on the tableau
    /// </summary>
    /// <param name="tableau">The simplex tableau</param>
    /// <param name="pivotCol">The pivot column index</param>
    /// <param name="pivotRow">The pivot row index</param>
    private void Pivot(double[,] tableau, int pivotCol, int pivotRow)
    {
        int rowCount = tableau.GetLength(0);
        int colCount = tableau.GetLength(1);
        
        // Store the pivot element
        double pivotElement = tableau[pivotRow, pivotCol];
        
        // Divide the pivot row by the pivot element
        for (int j = 0; j < colCount; j++)
        {
            tableau[pivotRow, j] /= pivotElement;
        }
        
        // Update other rows
        for (int i = 0; i < rowCount; i++)
        {
            if (i != pivotRow && Math.Abs(tableau[i, pivotCol]) > Epsilon)
            {
                double ratio = tableau[i, pivotCol];
                for (int j = 0; j < colCount; j++)
                {
                    tableau[i, j] -= ratio * tableau[pivotRow, j];
                }
            }
        }
    }
} // End of Pivot method
} // Closing class PrimalSimplexSolver
