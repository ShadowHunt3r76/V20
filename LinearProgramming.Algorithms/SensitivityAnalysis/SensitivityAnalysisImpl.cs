using System;
using System.Linq;
using LinearProgramming.Parsing;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;
using static LinearProgramming.Algorithms.SensitivityAnalysis.SensitivityOutputFormatter;

namespace LinearProgramming.Algorithms.SensitivityAnalysis
{
    public class SensitivityAnalysisImpl : BaseSensitivityAnalysis
    {
        private readonly double[,] _finalTableau;
        private readonly int[] _basis;
        private readonly Variable[] _variables;
        private bool _analysisPerformed = false;

        public SensitivityAnalysisImpl(
            double[,] finalTableau, 
            int[] basis, 
            Variable[] variables, 
            ConstraintType[] constraintTypes, 
            int numVariables, 
            int numConstraints)
            : base(
                variables?.Select(v => v.Name).ToArray() ?? Array.Empty<string>(),
                constraintTypes ?? Array.Empty<ConstraintType>(),
                numVariables,
                numConstraints)
        {
            _finalTableau = finalTableau ?? throw new ArgumentNullException(nameof(finalTableau));
            _basis = basis ?? throw new ArgumentNullException(nameof(basis));
            _variables = variables ?? throw new ArgumentNullException(nameof(variables));
        }

        /// <summary>
        /// Performs the sensitivity analysis on the linear programming solution.
        /// This implementation analyzes the final tableau to determine sensitivity information.
        /// </summary>
        public override void PerformAnalysis()
        {
            if (_analysisPerformed)
                return;

            // The analysis is performed on-demand when specific methods are called
            // (e.g., GetRHSRange, GetShadowPrice, etc.)
            // This is a no-op as the implementation is already set up for lazy evaluation
            
            _analysisPerformed = true;
        }

        /// <summary>
        /// Gets the range of values for a non-basic variable's coefficient in the objective function
        /// while maintaining the current optimal basis.
        /// </summary>
        public (double lowerBound, double upperBound) GetNonBasicVariableRange(int nonBasicVarIndex)
        {
            if (nonBasicVarIndex < 0 || nonBasicVarIndex >= _numVariables)
                throw new ArgumentOutOfRangeException(nameof(nonBasicVarIndex));

            if (_basis.Contains(nonBasicVarIndex))
                throw new ArgumentException("Variable is basic, not non-basic.", nameof(nonBasicVarIndex));

            int numRows = _finalTableau.GetLength(0);
            int numCols = _finalTableau.GetLength(1);
            int rhsCol = numCols - 1;

            double lowerBound = double.NegativeInfinity;
            double upperBound = double.PositiveInfinity;

            // For minimization problem, the reduced cost must remain non-negative to maintain optimality
            // For maximization, it should be non-positive (implementation depends on your problem type)
            // This example assumes a minimization problem
            for (int i = 1; i < numRows; i++)
            {
                double a_ij = _finalTableau[i, nonBasicVarIndex];
                if (Math.Abs(a_ij) > EPSILON)
                {
                    double ratio = _finalTableau[0, nonBasicVarIndex] / a_ij;
                    if (a_ij > 0 && ratio < upperBound)
                        upperBound = ratio;
                    else if (a_ij < 0 && ratio > lowerBound)
                        lowerBound = ratio;
                }
            }

            return (lowerBound, upperBound);
        }

        /// <summary>
        /// Gets the range of values for a basic variable's coefficient in the objective function
        /// while maintaining the current optimal basis.
        /// </summary>
        public (double lowerBound, double upperBound) GetBasicVariableRange(int basicVarIndex)
        {
            if (basicVarIndex < 0 || basicVarIndex >= _numVariables)
                throw new ArgumentOutOfRangeException(nameof(basicVarIndex));

            if (!_basis.Contains(basicVarIndex))
                throw new ArgumentException("Variable is non-basic, not basic.", nameof(basicVarIndex));

            int basisIndex = Array.IndexOf(_basis, basicVarIndex) + 1; // +1 for 1-based indexing in tableau
            int numCols = _finalTableau.GetLength(1);
            
            double lowerBound = double.NegativeInfinity;
            double upperBound = double.PositiveInfinity;

            for (int j = 0; j < numCols - 1; j++) // Exclude RHS column
            {
                if (_basis.Contains(j)) continue; // Skip basic variables
                
                double a_ij = _finalTableau[basisIndex, j];
                double c_j = _finalTableau[0, j];
                
                if (Math.Abs(a_ij) > EPSILON)
                {
                    double delta = -c_j / a_ij;
                    if (a_ij > 0 && delta < upperBound)
                        upperBound = delta;
                    else if (a_ij < 0 && delta > lowerBound)
                        lowerBound = delta;
                }
            }

            return (lowerBound, upperBound);
        }

        /// <summary>
        /// Gets the range of values for the right-hand side of a constraint
        /// while maintaining the current optimal basis.
        /// </summary>
        public override (double lowerBound, double upperBound) GetRHSRange(int constraintIndex)
        {
            if (constraintIndex < 0 || constraintIndex >= _numConstraints)
                throw new ArgumentOutOfRangeException(nameof(constraintIndex));

            int basisIndex = -1;
            for (int i = 0; i < _basis.Length; i++)
            {
                if (_basis[i] >= _numVariables && _basis[i] - _numVariables == constraintIndex)
                {
                    basisIndex = i + 1; // +1 for 1-based indexing in tableau
                    break;
                }
            }

            if (basisIndex == -1)
                throw new ArgumentException("Constraint not found in basis.", nameof(constraintIndex));

            int numCols = _finalTableau.GetLength(1);
            double currentRHS = _finalTableau[basisIndex, numCols - 1];
            
            double lowerBound = double.NegativeInfinity;
            double upperBound = double.PositiveInfinity;

            for (int j = 0; j < numCols - 1; j++) // Exclude RHS column
            {
                if (_basis.Contains(j)) continue; // Skip basic variables
                
                double a_ij = _finalTableau[basisIndex, j];
                double c_j = _finalTableau[0, j];
                
                if (Math.Abs(a_ij) > EPSILON)
                {
                    double delta = -c_j / a_ij;
                    if (a_ij > 0 && delta < upperBound)
                        upperBound = delta;
                    else if (a_ij < 0 && delta > lowerBound)
                        lowerBound = delta;
                }
            }

            // Adjust bounds based on constraint type
            switch (_constraintTypes[constraintIndex])
            {
                case ConstraintType.LessThanOrEqual:
                    return (currentRHS + lowerBound, currentRHS + upperBound);
                case ConstraintType.GreaterThanOrEqual:
                    return (currentRHS - upperBound, currentRHS - lowerBound);
                case ConstraintType.Equal:
                    // For equality constraints, the range is typically very narrow
                    return (currentRHS + lowerBound, currentRHS + upperBound);
                default:
                    throw new NotSupportedException($"Unsupported constraint type: {_constraintTypes[constraintIndex]}");
            }
        }

        /// <summary>
        /// Gets the shadow price (dual value) for a constraint.
        /// </summary>
        public override double GetShadowPrice(int constraintIndex)
        {
            if (constraintIndex < 0 || constraintIndex >= _numConstraints)
                throw new ArgumentOutOfRangeException(nameof(constraintIndex));

            int slackVarIndex = _numVariables + constraintIndex;
            return -_finalTableau[0, slackVarIndex]; // Negative of reduced cost for slack variable
        }

        /// <summary>
        /// Adds a new constraint to the optimal solution and determines if it affects the solution.
        /// </summary>
        public (bool isFeasible, double[] newSolution, double newObjectiveValue) 
            AddNewConstraint(double[] coefficients, double rhs, ConstraintType constraintType)
        {
            if (coefficients == null || coefficients.Length != _numVariables)
                throw new ArgumentException("Invalid number of coefficients.", nameof(coefficients));

            int numRows = _finalTableau.GetLength(0);
            int numCols = _finalTableau.GetLength(1);
            int rhsCol = numCols - 1;
            
            // Create a new tableau with an additional row
            var newTableau = new double[numRows + 1, numCols];
            
            // Copy existing tableau
            for (int i = 0; i < numRows; i++)
                for (int j = 0; j < numCols; j++)
                    newTableau[i, j] = _finalTableau[i, j];
            
            // Add new constraint row
            for (int j = 0; j < _numVariables; j++)
            {
                double sum = 0;
                for (int k = 0; k < _numVariables; k++)
                {
                    if (_basis.Contains(k))
                    {
                        int basisIndex = Array.IndexOf(_basis, k);
                        sum += coefficients[k] * _finalTableau[basisIndex + 1, j];
                    }
                }
                newTableau[numRows, j] = sum;
            }
            
            // Set the RHS for the new constraint
            double newRhs = rhs;
            for (int i = 0; i < _numVariables; i++)
            {
                if (_basis.Contains(i))
                {
                    int basisIndex = Array.IndexOf(_basis, i);
                    newRhs -= coefficients[i] * _finalTableau[basisIndex + 1, rhsCol];
                }
            }
            newTableau[numRows, rhsCol] = newRhs;
            
            // Check if the new constraint is violated
            bool isViolated = (constraintType == ConstraintType.LessThanOrEqual && newRhs < -EPSILON) ||
                             (constraintType == ConstraintType.GreaterThanOrEqual && newRhs > EPSILON) ||
                             (constraintType == ConstraintType.Equal && Math.Abs(newRhs) > EPSILON);
            
            if (!isViolated)
            {
                // Current solution is still feasible
                double[] currentSolution = new double[_numVariables];
                for (int i = 0; i < _numVariables; i++)
                {
                    if (_basis.Contains(i))
                    {
                        int basisIndex = Array.IndexOf(_basis, i);
                        currentSolution[i] = _finalTableau[basisIndex + 1, rhsCol];
                    }
                    else
                    {
                        currentSolution[i] = 0;
                    }
                }
                
                return (true, currentSolution, _finalTableau[0, rhsCol]);
            }
            
            // If we get here, we need to perform dual simplex iterations
            // This is a simplified version - a full implementation would handle all edge cases
            
            // Add slack/surplus variable
            int newSlackCol = numCols - 1; // Position for the new slack variable
            
            // Extend the tableau to include the new slack variable
            var extendedTableau = new double[numRows + 1, numCols + 1];
            for (int i = 0; i < numRows + 1; i++)
            {
                for (int j = 0; j < numCols; j++)
                {
                    if (j < newSlackCol)
                        extendedTableau[i, j] = newTableau[i, j];
                    else if (j > newSlackCol)
                        extendedTableau[i, j + 1] = newTableau[i, j];
                }
                
                // Add slack/surplus coefficient
                if (i == numRows) // New constraint row
                {
                    extendedTableau[i, newSlackCol] = (constraintType == ConstraintType.LessThanOrEqual) ? 1 : -1;
                }
                else if (i > 0) // Other constraint rows
                {
                    extendedTableau[i, newSlackCol] = 0;
                }
            }
            
            // Update the basis to include the new slack variable
            int newBasisVar = _numVariables + _numConstraints; // Index of the new slack variable
            var newBasis = _basis.ToList();
            newBasis.Add(newBasisVar);
            
            // Perform dual simplex iterations to restore feasibility
            // This is a simplified version - a full implementation would handle all edge cases
            
            // For now, return that the solution became infeasible
            return (false, null, double.NaN);
        }

        /// <summary>
        /// Adds a new variable (activity) to the optimal solution.
        /// </summary>
        public (bool isOptimal, double[] newSolution, double newObjectiveValue) 
            AddNewVariable(double objCoefficient, double[] constraintCoefficients)
        {
            if (constraintCoefficients == null || constraintCoefficients.Length != _numConstraints)
                throw new ArgumentException("Invalid number of constraint coefficients.", nameof(constraintCoefficients));

            int numRows = _finalTableau.GetLength(0);
            int numCols = _finalTableau.GetLength(1);
            int rhsCol = numCols - 1;
            
            // Calculate reduced cost for the new variable
            double reducedCost = objCoefficient;
            for (int i = 0; i < _numConstraints; i++)
            {
                if (_basis.Contains(_numVariables + i))
                {
                    int basisIndex = Array.IndexOf(_basis, _numVariables + i);
                    reducedCost -= constraintCoefficients[i] * _finalTableau[basisIndex + 1, 0];
                }
            }
            
            // If reduced cost is non-negative (for maximization) or non-positive (for minimization),
            // the current solution remains optimal
            // This assumes we're doing maximization - adjust the condition for minimization
            bool isOptimal = (reducedCost <= 0);
            
            if (isOptimal)
            {
                // Current solution remains optimal with the new variable at zero
                double[] currentSolution = new double[_numVariables + 1];
                for (int i = 0; i < _numVariables; i++)
                {
                    if (_basis.Contains(i))
                    {
                        int basisIndex = Array.IndexOf(_basis, i);
                        currentSolution[i] = _finalTableau[basisIndex + 1, rhsCol];
                    }
                    else
                    {
                        currentSolution[i] = 0;
                    }
                }
                currentSolution[_numVariables] = 0; // New variable at zero
                
                return (true, currentSolution, _finalTableau[0, rhsCol]);
            }
            
            // If we get here, we need to perform primal simplex iterations
            // to find the new optimal solution with the new variable
            
            // This is a simplified version - a full implementation would:
            // 1. Add a new column to the tableau
            // 2. Update the basis if needed
            // 3. Perform primal simplex iterations if necessary
            
            // For now, return that the solution is not optimal
            return (false, null, double.NaN);
        }
        
        /// <summary>
        /// Performs duality analysis on the linear programming problem.
        /// </summary>
        public DualityAnalysisResult AnalyzeDuality()
        {
            int numRows = _finalTableau.GetLength(0);
            int numCols = _finalTableau.GetLength(1);
            int rhsCol = numCols - 1;
            
            // Extract dual variables (shadow prices)
            double[] dualVariables = new double[_numConstraints];
            for (int i = 0; i < _numConstraints; i++)
            {
                int slackVarIndex = _numVariables + i;
                dualVariables[i] = -_finalTableau[0, slackVarIndex];
            }
            
            // Calculate dual objective value
            double dualObjective = 0;
            for (int i = 0; i < _numConstraints; i++)
            {
                dualObjective += dualVariables[i] * _finalTableau[i + 1, rhsCol];
            }
            
            // Get primal objective value
            double primalObjective = _finalTableau[0, rhsCol];
            
            // Check for strong/weak duality
            bool hasStrongDuality = Math.Abs(primalObjective - dualObjective) < EPSILON;
            
            return new DualityAnalysisResult
            {
                DualVariables = dualVariables,
                PrimalObjective = primalObjective,
                DualObjective = dualObjective,
                HasStrongDuality = hasStrongDuality,
                DualityGap = Math.Abs(primalObjective - dualObjective)
            };
        }
    }
    
    /// <summary>
    /// Represents the result of a duality analysis.
    /// </summary>
    public class DualityAnalysisResult
    {
        /// <summary>
        /// Gets the dual variables (shadow prices) for the constraints.
        /// </summary>
        public double[] DualVariables { get; set; }
        
        /// <summary>
        /// Gets the primal objective value.
        /// </summary>
        public double PrimalObjective { get; set; }
        
        /// <summary>
        /// Gets the dual objective value.
        /// </summary>
        public double DualObjective { get; set; }
        
        /// <summary>
        /// Gets a value indicating whether strong duality holds.
        /// </summary>
        public bool HasStrongDuality { get; set; }
        
        /// <summary>
        /// Gets the duality gap (absolute difference between primal and dual objectives).
        /// </summary>
        public double DualityGap { get; set; }
    }
    }
