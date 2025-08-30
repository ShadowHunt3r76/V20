using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.PrimalSimplex;
using static LinearProgramming.Parsing.ParsedLinearProgrammingModel;

// Use the outer CanonicalLinearProgrammingModel from LinearProgramming.Parsing
using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Implements the Revised Primal Simplex algorithm for solving linear programs.
    /// This implementation includes proper matrix operations, basis management, and Phase I/II method.
    /// </summary>
    public class RevisedPrimalSimplexSolver
    {
        // Event for iteration logging (commented out as it's not currently used)
        // public event Action<string> OnIteration;
        private const double EPSILON = MatrixUtils.Epsilon;
        private const double BIG_M = 1e6;
        private bool _useProductForm = true;
        private List<RevisedSimplexIteration> _iterations = new List<RevisedSimplexIteration>();
        private int _iterationCount = 0;
        private string[] _variableNames;
        private int _originalVarCount;
        private int _slackCount;
        private int _artificialCount;
        private int _totalVars;
        private double[,] _basisInverse;
        private double[] _simplexMultipliers;
        private double[] _reducedCosts;
        private double[] _currentSolution;
        private List<int> _basisIndices;
        private List<int> _nonBasisIndices;
        private double[][] _workingA;
        private double[] _workingB;
        private double[] _workingC;
        private double _currentObjective;
        private int _constraintCount;
        private bool _useBlandsRule = false; // Toggle between Dantzig's and Bland's rule
        
        #region Matrix Operations (Delegated to MatrixUtils)
        
        private double[,] MatrixInverse(double[,] matrix) => MatrixUtils.InvertMatrix(matrix, EPSILON);
        
        private double[] MatrixVectorMultiply(double[,] matrix, double[] vector) => 
            MatrixUtils.MatrixVectorMultiply(matrix, vector);
            
        private double[] VectorMatrixMultiply(double[] vector, double[,] matrix) => 
            MatrixUtils.VectorMatrixMultiply(vector, matrix);
            
        private double[,] ExtractBasisMatrix(double[][] A, List<int> basisIndices) => 
            MatrixUtils.ExtractBasisMatrix(A, basisIndices.ToArray(), basisIndices.Count);
            
        private double[] GetColumn(double[][] matrix, int colIndex) => 
            MatrixUtils.GetColumn(matrix, colIndex);
            
        #endregion

        #region Core Solver Methods
        
        private LinearProgramSolution InitializeSolution(CanonicalLinearProgrammingModel model)
        {
            _constraintCount = model.CoefficientMatrix.Length;
            _originalVarCount = model.CoefficientMatrix[0].Length;
            _slackCount = _constraintCount;
            _artificialCount = 0;
            
            // Count artificial variables needed for >= and = constraints
            foreach (var constraintType in model.ConstraintTypes)
            {
                if (constraintType == ConstraintType.GreaterThanOrEqual || 
                    constraintType == ConstraintType.Equal)
                {
                    _artificialCount++;
                }
            }
            
            _totalVars = _originalVarCount + _slackCount + _artificialCount;
            
            // Initialize variable names
            InitializeVariableNames();
            
            // Initialize working matrices and vectors
            InitializeWorkingMatrices(model);
            
            // Create initial basis
            InitializeBasis();
            
            // Create and return solution object
            var solution = new LinearProgramSolution
            {
                VariableNames = _variableNames,
                Status = "Initialized",
                UsingProductForm = _useProductForm
            };
            
            return solution;
        }
        
        private void InitializeVariableNames()
        {
            var names = new List<string>();
            
            // Original variables
            for (int i = 0; i < _originalVarCount; i++)
                names.Add($"x{i + 1}");
                
            // Slack variables
            for (int i = 0; i < _slackCount; i++)
                names.Add($"s{i + 1}");
                
            // Artificial variables
            for (int i = 0; i < _artificialCount; i++)
                names.Add($"a{i + 1}");
                
            _variableNames = names.ToArray();
        }
        
        private void InitializeWorkingMatrices(CanonicalLinearProgrammingModel model)
        {
            int m = _constraintCount;
            int n = _originalVarCount;
            
            // Initialize working A matrix (constraints)
            _workingA = new double[m][];
            for (int i = 0; i < m; i++)
            {
                _workingA[i] = new double[_totalVars];
                
                // Copy original coefficients
                for (int j = 0; j < n; j++)
                {
                    _workingA[i][j] = model.CoefficientMatrix[i][j];
                }
                
                // Add slack/surplus variables
                int slackPos = _originalVarCount + i;
                if (model.ConstraintTypes[i] == ConstraintType.LessThanOrEqual)
                {
                    _workingA[i][slackPos] = 1.0; // Slack variable
                }
                else if (model.ConstraintTypes[i] == ConstraintType.GreaterThanOrEqual)
                {
                    _workingA[i][slackPos] = -1.0; // Surplus variable
                }
                // For = constraints, no slack/surplus variable is added
            }
            
            // Initialize working b vector (RHS)
            _workingB = new double[m];
            Array.Copy(model.RightHandSide, _workingB, m);
            
            // Initialize working c vector (objective)
            _workingC = new double[_totalVars];
            Array.Copy(model.ObjectiveCoefficients, _workingC, n);
            
            // Set coefficients for artificial variables in the objective (for Phase I)
            for (int i = _originalVarCount + _slackCount; i < _totalVars; i++)
            {
                _workingC[i] = BIG_M; // Big M method for artificial variables
            }
        }
        
        private void InitializeBasis()
        {
            _basisIndices = new List<int>();
            _nonBasisIndices = new List<int>();
            
            // Add slack variables to basis first
            for (int i = 0; i < _constraintCount; i++)
            {
                int slackPos = _originalVarCount + i;
                _basisIndices.Add(slackPos);
            }
            
            // Add artificial variables to basis where needed
            int artVar = _originalVarCount + _slackCount;
            for (int i = 0; i < _artificialCount; i++)
            {
                _basisIndices.Add(artVar + i);
            }
            
            // All other variables are non-basic
            for (int i = 0; i < _totalVars; i++)
            {
                if (!_basisIndices.Contains(i))
                {
                    _nonBasisIndices.Add(i);
                }
            }
            
            // Initialize basis inverse
            var basisMatrix = ExtractBasisMatrix(_workingA, _basisIndices);
            _basisInverse = MatrixInverse(basisMatrix);
            
            if (_basisInverse == null)
            {
                throw new InvalidOperationException("Initial basis matrix is singular");
            }
            
            // Initialize current solution
            UpdateSolution();
        }
        
        private (int entering, double[] reducedCosts) PriceOut()
        {
            int m = _constraintCount;
            int n = _totalVars;
            
            // Compute simplex multipliers: π = cB^T * B^(-1)
            double[] cB = new double[m];
            for (int i = 0; i < m; i++)
            {
                cB[i] = _workingC[_basisIndices[i]];
            }
            _simplexMultipliers = VectorMatrixMultiply(cB, _basisInverse);
            
            // Compute reduced costs for non-basic variables
            _reducedCosts = new double[n];
            int entering = -1;
            double mostNegative = 0.0;
            
            foreach (int j in _nonBasisIndices)
            {
                double[] Aj = GetColumn(_workingA, j);
                _reducedCosts[j] = _workingC[j] - DotProduct(_simplexMultipliers, Aj);
                
                // Check for most negative reduced cost (Dantzig's rule)
                if (_reducedCosts[j] < mostNegative - EPSILON)
                {
                    mostNegative = _reducedCosts[j];
                    entering = j;
                }
                else if (_useBlandsRule && 
                         Math.Abs(_reducedCosts[j] - mostNegative) < EPSILON && 
                         j < entering)
                {
                    // Bland's rule: choose variable with smallest index in case of tie
                    entering = j;
                }
            }
            
            return (entering, _reducedCosts);
        }
        
        private (int leaving, double minRatio) RatioTest(int entering)
        {
            int m = _constraintCount;
            double[] Aj = GetColumn(_workingA, entering);
            double[] d = MatrixVectorMultiply(_basisInverse, Aj);
            
            int leaving = -1;
            double minRatio = double.MaxValue;
            
            for (int i = 0; i < m; i++)
            {
                if (d[i] > EPSILON) // Only consider positive denominators
                {
                    double ratio = _currentSolution[i] / d[i];
                    
                    if (ratio < minRatio - EPSILON || 
                        (Math.Abs(ratio - minRatio) < EPSILON && _basisIndices[i] < _basisIndices[leaving]))
                    {
                        minRatio = ratio;
                        leaving = i;
                    }
                }
            }
            
            return (leaving, minRatio);
        }
        
        private void UpdateBasis(int entering, int leaving)
        {
            // Update basis and non-basis indices
            int leavingVar = _basisIndices[leaving];
            _basisIndices[leaving] = entering;
            _nonBasisIndices.Remove(entering);
            _nonBasisIndices.Add(leavingVar);
            
            if (_useProductForm)
            {
                UpdateBasisProductForm(entering, leaving);
            }
            else
            {
                UpdateBasisPriceOut(entering, leaving);
            }
        }
        
        private void UpdateBasisProductForm(int entering, int leaving)
        {
            // Validate inputs
            if (entering < 0 || entering >= _totalVars)
                throw new ArgumentOutOfRangeException(nameof(entering), "Entering variable index is out of range");
                
            if (leaving < 0 || leaving >= _constraintCount)
                throw new ArgumentOutOfRangeException(nameof(leaving), "Leaving variable index is out of range");
                
            if (_basisInverse == null)
                throw new InvalidOperationException("Basis inverse is not initialized");
                
            if (_workingA == null)
                throw new InvalidOperationException("Working matrix A is not initialized");
            // Get the entering column in the original basis
            double[] Aj = GetColumn(_workingA, entering);
            double[] d = MatrixVectorMultiply(_basisInverse, Aj);
            
            // Create eta matrix for the update
            int m = _constraintCount;
            double[,] eta = new double[m, m];
            
            // Initialize as identity matrix
            for (int i = 0; i < m; i++)
                eta[i, i] = 1.0;
                
            // Update the eta matrix
            for (int i = 0; i < m; i++)
            {
                if (i == leaving)
                    eta[i, leaving] = 1.0 / d[leaving];
                else
                    eta[i, leaving] = -d[i] / d[leaving];
            }
            
            // Update basis inverse: B_new^(-1) = E * B_old^(-1)
            _basisInverse = MatrixMultiply(eta, _basisInverse);
            
            DisplayProductFormUpdate(entering, leaving, d);
        }
        
        private void DisplayProductFormUpdate(int entering, int leaving, double[] etaVector)
        {
            Console.WriteLine("\n=== Product Form Update ===");
            Console.WriteLine($"Entering: {_variableNames[entering]}, Leaving: {_basisIndices[leaving]} ({_variableNames[_basisIndices[leaving]]})");
            
            // Display eta vector
            Console.WriteLine("\nEta Vector (B^-1 * A_entering):");
            for (int i = 0; i < etaVector.Length; i++)
            {
                if (Math.Abs(etaVector[i]) > EPSILON)
                {
                    Console.WriteLine($"  eta[{i+1}] = {etaVector[i]:F6}");
                }
            }
            
            // Display eta matrix
            Console.WriteLine("\nEta Matrix (E):");
            double[,] etaMatrix = CreateIdentityMatrix(_constraintCount);
            for (int i = 0; i < _constraintCount; i++)
            {
                if (i == leaving)
                    etaMatrix[leaving, leaving] = 1.0 / etaVector[leaving];
                else
                    etaMatrix[i, leaving] = -etaVector[i] / etaVector[leaving];
                
                // Display row
                Console.Write("  [");
                for (int j = 0; j < _constraintCount; j++)
                {
                    if (j == leaving) // Highlight the column being updated
                        Console.ForegroundColor = ConsoleColor.Yellow;
                        
                    Console.Write($"{etaMatrix[i, j],8:F4} ");
                    
                    if (j == leaving)
                        Console.ResetColor();
                }
                Console.WriteLine("]");
            }
            Console.ResetColor();
        }
        
        private void UpdateBasisPriceOut(int entering, int leaving)
        {
            // Validate inputs
            if (entering < 0 || entering >= _totalVars)
                throw new ArgumentOutOfRangeException(nameof(entering), "Entering variable index is out of range");
                
            if (leaving < 0 || leaving >= _constraintCount)
                throw new ArgumentOutOfRangeException(nameof(leaving), "Leaving variable index is out of range");
                
            if (_basisIndices == null || _basisIndices.Count == 0)
                throw new InvalidOperationException("Basis indices are not initialized");
                
            if (_workingA == null)
                throw new InvalidOperationException("Working matrix A is not initialized");
            // For Price Out method, we'll recompute the full inverse
            // This is less efficient but more numerically stable
            var basisMatrix = ExtractBasisMatrix(_workingA, _basisIndices);
            _basisInverse = MatrixInverse(basisMatrix);
            
            DisplayPriceOutIteration(_iterationCount, entering, leaving, _reducedCosts);
        }
        
        private void DisplayPriceOutIteration(int iter, int entering, int leaving, double[] reducedCosts)
        {
            Console.WriteLine("\n=== Price Out Iteration ===");
            Console.WriteLine($"Iteration: {iter}, Entering: {_variableNames[entering]}");
            
            if (leaving >= 0)
                Console.WriteLine($"Leaving: {_basisIndices[leaving]} ({_variableNames[_basisIndices[leaving]]})");
                
            // Display reduced costs
            Console.WriteLine("\nReduced Costs (c_j - y*A_j):");
            for (int j = 0; j < reducedCosts.Length; j++)
            {
                if (!_basisIndices.Contains(j) && Math.Abs(reducedCosts[j]) > EPSILON)
                {
                    string varName = j < _variableNames.Length ? _variableNames[j] : $"s{j - _originalVarCount + 1}";
                    Console.WriteLine($"  {varName}: {reducedCosts[j]:F6}");
                }
            }
            
            // Display current basis and solution
            Console.WriteLine("\nCurrent Basis:");
            for (int i = 0; i < _basisIndices.Count; i++)
            {
                string varName = _basisIndices[i] < _variableNames.Length ? 
                    _variableNames[_basisIndices[i]] : $"s{_basisIndices[i] - _originalVarCount + 1}";
                Console.WriteLine($"  {varName} = {_currentSolution[_basisIndices[i]]:F4}");
            }
            
            Console.WriteLine($"\nCurrent Objective: {_currentObjective:F4}\n");
        }
        
        private double[,] MatrixMultiply(double[,] a, double[,] b)
        {
            int m = a.GetLength(0);
            int n = b.GetLength(1);
            int p = a.GetLength(1);
            
            double[,] result = new double[m, n];
            
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < p; k++)
                    {
                        sum += a[i, k] * b[k, j];
                    }
                    result[i, j] = sum;
                }
            }
            
            return result;
        }
        
        private void UpdateSolution()
        {
            try
            {
                if (_basisInverse == null)
                    throw new InvalidOperationException("Basis inverse is not initialized");
                    
                if (_workingB == null)
                    throw new InvalidOperationException("RHS vector is not initialized");
                    
                // Check for numerical stability
                for (int i = 0; i < _constraintCount; i++)
                {
                    if (double.IsNaN(_workingB[i]) || double.IsInfinity(_workingB[i]))
                        throw new InvalidOperationException("Numerical instability detected in RHS values");
                }
                
                // Solve B * xB = b
                _currentSolution = MatrixVectorMultiply(_basisInverse, _workingB);
                
                // Update objective value
                double[] cB = new double[_constraintCount];
                for (int i = 0; i < _constraintCount; i++)
                {
                    cB[i] = _workingC[_basisIndices[i]];
                }
                
                _currentObjective = DotProduct(cB, _currentSolution);
            }
            catch (Exception ex)
            {
                throw new InvalidOperationException("Error updating solution: " + ex.Message, ex);
            }
        }
        
        private RevisedSimplexIteration CreateIteration()
        {
            var iteration = new RevisedSimplexIteration
            {
                Iteration = _iterationCount,
                Type = _useProductForm ? "ProductForm" : "PriceOut",
                BasicSolution = (double[])_currentSolution.Clone(),
                SimplexMultipliers = (double[])_simplexMultipliers?.Clone(),
                ReducedCosts = (double[])_reducedCosts?.Clone(),
                ObjectiveValue = _currentObjective,
                Status = "In Progress"
            };
            
            // Set basis indices through the BasisIndices property
            iteration.BasisIndices = new List<int>(_basisIndices);
            
            return iteration;
        }
        
        #endregion

        #region Phase I and II Methods
        
        private (bool IsFeasible, List<int> BasisIndices, List<int> NonBasisIndices, double[,] BasisInverse, double[] Solution, double ObjectiveValue) 
            RunPhaseI(CanonicalLinearProgrammingModel model)
        {
            _iterationCount = 0;
            bool isOptimal = false;
            
            while (!isOptimal && _iterationCount < 1000) // Safety limit
            {
                _iterationCount++;
                
                // Create iteration record
                var iteration = CreateIteration();
                
                // Pricing out: find entering variable
                var (entering, reducedCosts) = PriceOut();
                
                // Check for optimality (all reduced costs >= 0)
                if (entering == -1 || reducedCosts[entering] >= -EPSILON)
                {
                    isOptimal = true;
                    
                    // Check if we have artificial variables in the basis with non-zero values
                    bool hasArtificialInBasis = false;
                    for (int i = 0; i < _constraintCount; i++)
                    {
                        if (_basisIndices[i] >= _originalVarCount + _slackCount && 
                            Math.Abs(_currentSolution[i]) > EPSILON)
                        {
                            hasArtificialInBasis = true;
                            break;
                        }
                    }
                    
                    if (hasArtificialInBasis)
                    {
                        iteration.Status = "Infeasible (artificial variables in basis)";
                        _iterations.Add(iteration);
                        return (false, null, null, null, null, 0);
                    }
                    
                    iteration.Status = "Phase I Complete";
                    _iterations.Add(iteration);
                    break;
                }
                
                // Ratio test to find leaving variable
                var (leaving, minRatio) = RatioTest(entering);
                
                if (leaving == -1)
                {
                    iteration.Status = "Unbounded in Phase I";
                    _iterations.Add(iteration);
                    return (false, null, null, null, null, 0);
                }
                
                // Update iteration with pivot information
                iteration.EnteringVariable = entering;
                iteration.LeavingVariable = _basisIndices[leaving];
                iteration.PivotElement = GetColumn(_workingA, entering)[leaving];
                iteration.Status = "Pivoting";
                
                // Update basis and solution
                UpdateBasis(entering, leaving);
                UpdateSolution();
                
                // Add completed iteration
                _iterations.Add(iteration);
            }
            
            return (true, _basisIndices, _nonBasisIndices, _basisInverse, _currentSolution, _currentObjective);
        }
        
        private (bool IsOptimal, double ObjectiveValue) RunPhaseII(CanonicalLinearProgrammingModel model)
        {
            _iterationCount = 0;
            bool isOptimal = false;
            
            while (!isOptimal && _iterationCount < 1000) // Safety limit
            {
                _iterationCount++;
                
                // Create iteration record
                var iteration = CreateIteration();
                
                // Pricing out: find entering variable
                var (entering, reducedCosts) = PriceOut();
                
                // Check for optimality (all reduced costs >= 0 for maximization)
                if (entering == -1 || reducedCosts[entering] >= -EPSILON)
                {
                    isOptimal = true;
                    iteration.Status = "Optimal Solution Found";
                    _iterations.Add(iteration);
                    break;
                }
                
                // Ratio test to find leaving variable
                var (leaving, minRatio) = RatioTest(entering);
                
                if (leaving == -1)
                {
                    iteration.Status = "Unbounded";
                    _iterations.Add(iteration);
                    return (false, double.PositiveInfinity);
                }
                
                // Update iteration with pivot information
                iteration.EnteringVariable = entering;
                iteration.LeavingVariable = _basisIndices[leaving];
                iteration.PivotElement = GetColumn(_workingA, entering)[leaving];
                iteration.Status = "Pivoting";
                
                // Update basis and solution
                UpdateBasis(entering, leaving);
                UpdateSolution();
                
                // Add completed iteration
                _iterations.Add(iteration);
            }
            
            return (isOptimal, _currentObjective);
        }
        
        private void RemoveArtificialVariables()
        {
            // Try to remove artificial variables from the basis by replacing them with
            // non-basic original or slack variables if possible
            for (int i = 0; i < _constraintCount; i++)
            {
                int basisVar = _basisIndices[i];
                
                // Skip if not an artificial variable
                if (basisVar < _originalVarCount + _slackCount)
                    continue;
                    
                // Look for a non-basic original or slack variable to enter the basis
                foreach (int nonBasicVar in _nonBasisIndices)
                {
                    if (nonBasicVar < _originalVarCount + _slackCount) // Original or slack variable
                    {
                        double[] Aj = GetColumn(_workingA, nonBasicVar);
                        double[] d = MatrixVectorMultiply(_basisInverse, Aj);
                        
                        // If this variable can replace the artificial variable
                        if (Math.Abs(d[i]) > EPSILON)
                        {
                            // Update basis
                            _basisIndices[i] = nonBasicVar;
                            _nonBasisIndices.Remove(nonBasicVar);
                            _nonBasisIndices.Add(basisVar);
                            
                            // Update basis inverse using the product form
                            UpdateBasisProductForm(nonBasicVar, i);
                            break;
                        }
                    }
                }
            }
            
            // Update solution after basis changes
            UpdateSolution();
        }
        
        #endregion

        public LinearProgramSolution Solve(CanonicalLinearProgrammingModel model)
        {
            // Input validation
            if (model == null)
                throw new ArgumentNullException(nameof(model), "Model cannot be null");
                
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
            
            // Check RHS values
            foreach (var rhs in model.RightHandSide)
            {
                if (double.IsNaN(rhs) || double.IsInfinity(rhs))
                    throw new ArgumentException("RHS values must be finite numbers");
            }
            // Initialize solution
            var solution = new LinearProgramSolution
            {
                Status = "In Progress",
                VariableNames = model.VariableNames?.ToArray() ?? 
                    Enumerable.Range(1, model.ObjectiveCoefficients.Length).Select(i => $"x{i}").ToArray()
            };
            
            // Display problem information
            Console.WriteLine("\n=== Solving Linear Program with Revised Primal Simplex ===");
            DisplayCanonicalForm(model);
            Console.WriteLine("\n=== Starting Revised Primal Simplex Method ===");
            
            // Initialize solution and problem data
            var currentSolution = InitializeSolution(model);
            bool phaseIneeded = _artificialCount > 0;
            
            // Setup iteration logging after currentSolution is initialized
            // Commented out as we're not using the OnIteration event anymore
            // OnIteration += (message) => 
            // {
            //     Console.WriteLine(message);
            //     // currentSolution.Iterations.Add(new TableauIteration 
            //     // { 
            //     //     Iteration = _iterationCount,
            //     //     Description = message
            //     // });
            // };
            
            try
            {
                
                // Phase I: Find initial feasible solution if needed
                if (phaseIneeded)
                {
                    currentSolution.Status = "Running Phase I";
                    var phaseIResult = RunPhaseI(model);
                    
                    if (!phaseIResult.IsFeasible)
                    {
                        currentSolution.Status = "Infeasible";
                        currentSolution.IterationDetails = _iterations;
                        return currentSolution;
                    }
                    
                    // Update basis and solution for Phase II
                    _basisIndices = phaseIResult.BasisIndices;
                    _nonBasisIndices = phaseIResult.NonBasisIndices;
                    _basisInverse = phaseIResult.BasisInverse;
                    _currentSolution = phaseIResult.Solution;
                    _currentObjective = phaseIResult.ObjectiveValue;
                    
                    // Remove artificial variables from basis if possible
                    RemoveArtificialVariables();
                    
                    // Switch to original objective for Phase II
                    InitializeWorkingMatrices(model);
                    currentSolution.Status = "Phase I Complete, Starting Phase II";
                }
                
                // Phase II: Optimize with original objective
                currentSolution.Status = "Running Phase II";
                var phaseIIResult = RunPhaseII(model);
                
                // Prepare final solution
                currentSolution.ObjectiveValue = _currentObjective;
                currentSolution.IterationDetails = _iterations;
                currentSolution.Status = phaseIIResult.IsOptimal ? "Optimal" : "Unbounded or Infeasible";
                
                // Extract solution values
                currentSolution.SolutionVector = new double[_originalVarCount];
                for (int i = 0; i < _constraintCount; i++)
                {
                    if (_basisIndices[i] < _originalVarCount) // Only original variables
                    {
                        currentSolution.SolutionVector[_basisIndices[i]] = _currentSolution[i];
                    }
                }
                
                // Store basis and reduced costs
                currentSolution.Basis = _basisIndices.ToArray();
                currentSolution.ReducedCosts = _reducedCosts;
                currentSolution.SimplexMultipliers = _simplexMultipliers;
                
                return currentSolution;
            }
            catch (Exception ex)
            {
                return new LinearProgramSolution
                {
                    Status = $"Error: {ex.Message}",
                    IterationDetails = _iterations
                };
            }
        }
        
        /// <summary>
        /// Finds initial basis and determines if Phase I is needed
        /// </summary>
        private (List<int> Basis, bool NeedsPhaseI, double[][] WorkingA, double[] WorkingC, double[] WorkingB) 
            FindInitialBasis(double[][] A, double[] b, VariableType[] varTypes, int m, int n)
        {
            var basis = new List<int>();
            bool needsPhaseI = false;
            var workingA = new double[m][];
            var workingC = new List<double>();
            var workingB = new double[m];
            
            // Copy original matrix and objective
            for (int i = 0; i < m; i++)
            {
                workingA[i] = new double[n + m]; // Allow space for artificial variables
                Array.Copy(A[i], workingA[i], n);
                workingB[i] = b[i];
            }
            
            // Try to find identity columns in the original matrix
            var identityColumns = FindIdentityColumns(A, m, n);
            
            if (identityColumns.Count == m)
            {
                // Found natural basis
                basis = identityColumns;
                needsPhaseI = false;
                
                // Extend workingA to match expected size
                for (int i = 0; i < m; i++)
                {
                    Array.Resize(ref workingA[i], n);
                }
            }
            else
            {
                // Need artificial variables
                needsPhaseI = true;
                int artificialIndex = n;
                
                // Add artificial variables where needed
                for (int i = 0; i < m; i++)
                {
                    if (i < identityColumns.Count)
                    {
                        basis.Add(identityColumns[i]);
                    }
                    else
                    {
                        // Add artificial variable
                        basis.Add(artificialIndex);
                        workingA[i][artificialIndex] = 1.0;
                        artificialIndex++;
                    }
                }
                
                // Adjust working matrix size
                int totalVars = artificialIndex;
                for (int i = 0; i < m; i++)
                {
                    Array.Resize(ref workingA[i], totalVars);
                }
                
                // Create Phase I objective (minimize sum of artificial variables)
                workingC = new double[totalVars].ToList();
                for (int j = n; j < totalVars; j++)
                {
                    workingC[j] = 1.0; // Minimize artificial variables
                }
            }
            
            return (basis, needsPhaseI, workingA, workingC.ToArray(), workingB);
        }
        
        /// <summary>
        /// Finds columns that form identity matrix
        /// </summary>
        private List<int> FindIdentityColumns(double[][] A, int m, int n)
        {
            var identityColumns = new List<int>();
            var usedRows = new bool[m];
            
            // Look for columns that have exactly one 1 and rest zeros
            for (int j = 0; j < n; j++)
            {
                int oneCount = 0;
                int oneRow = -1;
                bool isCandidate = true;
                
                for (int i = 0; i < m; i++)
                {
                    if (Math.Abs(A[i][j] - 1.0) < EPSILON)
                    {
                        oneCount++;
                        oneRow = i;
                    }
                    else if (Math.Abs(A[i][j]) > EPSILON)
                    {
                        isCandidate = false;
                        break;
                    }
                }
                
                if (isCandidate && oneCount == 1 && !usedRows[oneRow])
                {
                    identityColumns.Add(j);
                    usedRows[oneRow] = true;
                }
            }
            
            return identityColumns.OrderBy(x => x).ToList();
        }
        
        /// <summary>
        /// Solves Phase I to find initial feasible solution
        /// </summary>
        private (string Status, List<int> Basis, double[] BasicSolution) 
            SolvePhaseI(double[][] A, double[] b, List<int> basis, int m)
        {
            // Create Phase I objective: minimize sum of artificial variables
            int n = A[0].Length;
            double[] phaseIObjective = new double[n];
            
            // Set coefficients for artificial variables to 1
            for (int j = 0; j < basis.Count; j++)
            {
                if (basis[j] >= n - m) // Artificial variable
                {
                    phaseIObjective[basis[j]] = 1.0;
                }
            }
            
            // Solve Phase I problem
            var result = SolveSimplex(A, b, phaseIObjective, basis, m, true);
            
            if (result.Status == "Optimal" && Math.Abs(result.ObjectiveValue) < EPSILON)
            {
                return ("Optimal", result.Basis, result.BasicSolution);
            }
            else
            {
                return ("Infeasible", basis, null);
            }
        }
        
        /// <summary>
        /// Solves Phase II with original objective
        /// </summary>
        private LinearProgramSolution SolvePhaseII(double[][] A, double[] b, double[] c, List<int> basis, int m, int n)
        {
            var result = SolveSimplex(A, b, c, basis, m, false);
            
            var solution = new LinearProgramSolution
            {
                Status = result.Status,
                ObjectiveValue = result.ObjectiveValue
            };
            
            if (result.Status == "Optimal")
            {
                solution.SolutionVector = ConstructFullSolution(result.Basis, result.BasicSolution, n);
            }
            
            return solution;
        }
        
        /// <summary>
        /// Core simplex algorithm implementation
        /// </summary>
        private (string Status, List<int> Basis, double[] BasicSolution, double ObjectiveValue)
            SolveSimplex(double[][] A, double[] b, double[] c, List<int> basis, int m, bool isPhaseI)
        {
            int n = A[0].Length;
            int maxIter = 1000;
            int iter = 0;
            
            while (iter < maxIter)
            {
                // Extract basis matrix B and compute inverse
                double[,] B = ExtractBasisMatrix(A, basis, m);
                double[,] BInverse = InvertMatrix(B);
                
                if (BInverse == null)
                {
                    // Try to repair basis by finding alternative basic variable
                    var newBasis = TryRepairBasis(A, basis, m, n);
                    if (newBasis != null)
                    {
                        basis = newBasis;
                        continue;
                    }
                    return ("Infeasible", basis, null, 0.0);
                }
                
                // Step 1: Compute basic feasible solution x_B = B^(-1) * b
                double[] xB = MultiplyMatrixVector(BInverse, b);
                
                // Check feasibility with tolerance for numerical errors
                bool feasible = true;
                for (int i = 0; i < xB.Length; i++)
                {
                    if (xB[i] < -EPSILON)
                    {
                        feasible = false;
                        break;
                    }
                }
                
                if (!feasible)
                {
                    return ("Infeasible", basis, null, 0.0);
                }
                
                // Step 2: Compute simplex multipliers ? = c_B^T * B^(-1)
                double[] cB = new double[basis.Count];
                for (int i = 0; i < basis.Count; i++) 
                    cB[i] = c[basis[i]];
                
                double[] pi = MultiplyVectorMatrix(cB, BInverse);
                
                // Step 3: Compute reduced costs for non-basic variables
                int entering = -1;
                double mostNegative = isPhaseI ? EPSILON : 0; // Different tolerance for Phase I
                
                for (int j = 0; j < n; j++)
                {
                    if (!basis.Contains(j))
                    {
                        // Extract column j from A
                        double[] Aj = ExtractColumn(A, j, m);
                        // Compute reduced cost: c_j - ?^T * A_j
                        double reducedCost = c[j] - DotProduct(pi, Aj);
                        
                        if (reducedCost < mostNegative)
                        {
                            mostNegative = reducedCost;
                            entering = j;
                        }
                    }
                }
                
                // Step 4: Check optimality
                if (entering == -1)
                {
                    // Optimal solution found
                    double objectiveValue = DotProduct(cB, xB);
                    return ("Optimal", basis, xB, objectiveValue);
                }
                
                // Step 5: Compute direction vector d = B^(-1) * A_entering
                double[] AEntering = ExtractColumn(A, entering, m);
                double[] d = MultiplyMatrixVector(BInverse, AEntering);
                
                // Step 6: Ratio test to find leaving variable
                int leaving = -1;
                double minRatio = double.MaxValue;
                
                for (int i = 0; i < d.Length; i++)
                {
                    if (d[i] > EPSILON)
                    {
                        double ratio = xB[i] / d[i];
                        if (ratio >= 0 && ratio < minRatio)
                        {
                            minRatio = ratio;
                            leaving = i;
                        }
                    }
                }
                
                if (leaving == -1)
                {
                    return ("Unbounded", basis, xB, double.PositiveInfinity);
                }
                
                // Step 7: Update basis
                basis[leaving] = entering;
                
                iter++;
            }
            
            return ("Max iterations reached", basis, null, 0.0);
        }
        
        /// <summary>
        /// Attempts to repair a singular basis by finding alternative basic variables
        /// </summary>
        private List<int> TryRepairBasis(double[][] A, List<int> currentBasis, int m, int n)
        {
            // Try replacing each basis variable with non-basic variables
            for (int i = 0; i < currentBasis.Count; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (!currentBasis.Contains(j))
                    {
                        var testBasis = new List<int>(currentBasis);
                        testBasis[i] = j;
                        
                        var B = ExtractBasisMatrix(A, testBasis, m);
                        if (InvertMatrix(B) != null)
                        {
                            return testBasis;
                        }
                    }
                }
            }
            return null;
        }
        
        #region Matrix Operations
        
        /// <summary>
        /// Extracts the basis matrix B from the coefficient matrix A
        /// </summary>
        private double[,] ExtractBasisMatrix(double[][] A, List<int> basis, int m)
        {
            double[,] B = new double[m, m];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    B[i, j] = A[i][basis[j]];
                }
            }
            return B;
        }
        
        /// <summary>
        /// Extracts column j from matrix A
        /// </summary>
        private double[] ExtractColumn(double[][] A, int j, int m)
        {
            double[] column = new double[m];
            for (int i = 0; i < m; i++)
            {
                column[i] = A[i][j];
            }
            return column;
        }
        
        /// <summary>
        /// Matrix inversion using Gauss-Jordan elimination
        /// Returns null if matrix is singular
        /// </summary>
        private double[,] InvertMatrix(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            double[,] augmented = new double[n, 2 * n];
            
            // Create augmented matrix [A|I]
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    augmented[i, j] = matrix[i, j];
                    augmented[i, j + n] = (i == j) ? 1.0 : 0.0;
                }
            }
            
            // Gauss-Jordan elimination
            for (int i = 0; i < n; i++)
            {
                // Find pivot
                int pivotRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(augmented[k, i]) > Math.Abs(augmented[pivotRow, i]))
                        pivotRow = k;
                }
                
                // Check for singular matrix
                if (Math.Abs(augmented[pivotRow, i]) < EPSILON)
                    return null;
                
                // Swap rows if necessary
                if (pivotRow != i)
                {
                    for (int j = 0; j < 2 * n; j++)
                    {
                        double temp = augmented[i, j];
                        augmented[i, j] = augmented[pivotRow, j];
                        augmented[pivotRow, j] = temp;
                    }
                }
                
                // Make diagonal element 1
                double pivot = augmented[i, i];
                for (int j = 0; j < 2 * n; j++)
                {
                    augmented[i, j] /= pivot;
                }
                
                // Eliminate column
                for (int k = 0; k < n; k++)
                {
                    if (k != i)
                    {
                        double factor = augmented[k, i];
                        for (int j = 0; j < 2 * n; j++)
                        {
                            augmented[k, j] -= factor * augmented[i, j];
                        }
                    }
                }
            }
            
            // Extract inverse matrix
            double[,] inverse = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    inverse[i, j] = augmented[i, j + n];
                }
            }
            
            return inverse;
        }
        
        /// <summary>
        /// Matrix-vector multiplication: result = A * x
        /// </summary>
        private double[] MultiplyMatrixVector(double[,] A, double[] x)
        {
            int m = A.GetLength(0);
            int n = A.GetLength(1);
            double[] result = new double[m];
            
            for (int i = 0; i < m; i++)
            {
                result[i] = 0;
                for (int j = 0; j < n; j++)
                {
                    result[i] += A[i, j] * x[j];
                }
            }
            
            return result;
        }
        
        /// <summary>
        /// Vector-matrix multiplication: result = x^T * A
        /// </summary>
        private double[] MultiplyVectorMatrix(double[] x, double[,] A)
        {
            int m = A.GetLength(0);
            int n = A.GetLength(1);
            double[] result = new double[n];
            
            for (int j = 0; j < n; j++)
            {
                result[j] = 0;
                for (int i = 0; i < m; i++)
                {
                    result[j] += x[i] * A[i, j];
                }
            }
            
            return result;
        }
        
        // DotProduct is implemented in MatrixUtils
        
        /// <summary>
        /// Constructs the full solution vector from basic variables
        /// </summary>
        private double[] ConstructFullSolution(List<int> basis, double[] xB, int n)
        {
            double[] x = new double[n];
            
            // Initialize all variables to 0 (non-basic variables)
            for (int i = 0; i < n; i++)
                x[i] = 0;
            
            // Set basic variables (only for original variables, not artificial)
            for (int i = 0; i < basis.Count; i++)
            {
                if (basis[i] < n) // Only original variables
                {
                    x[basis[i]] = xB[i];
                }
            }
            
            return x;
        }
        
        #endregion
        
        #region Helper Methods
        
        private double DotProduct(double[] a, double[] b)
        {
            return MatrixUtils.DotProduct(a, b);
        }
        
        private double[,] CreateIdentityMatrix(int size)
        {
            return MatrixUtils.CreateIdentityMatrix(size);
        }
        
        private void DisplayCanonicalForm(CanonicalLinearProgrammingModel model)
        {
            Console.WriteLine("Canonical Form:");
            
            // Default to maximization if we can't determine the objective type
            // In the canonical form, we'll assume maximization unless we have information to the contrary
            // The actual optimization direction is handled by the solver
            Console.WriteLine("Objective: Maximize");
            
            // Display objective function
            if (model.ObjectiveCoefficients != null)
            {
                Console.Write("  ");
                for (int j = 0; j < model.ObjectiveCoefficients.Length; j++)
                {
                    if (j > 0 && model.ObjectiveCoefficients[j] >= 0)
                        Console.Write(" + ");
                    else if (model.ObjectiveCoefficients[j] < 0)
                        Console.Write(" - ");
                    
                    Console.Write($"{Math.Abs(model.ObjectiveCoefficients[j]):F2}x{j + 1}");
                }
                Console.WriteLine("\n");
            }
            
            // Display constraints
            Console.WriteLine("Subject to:");
            for (int i = 0; i < model.CoefficientMatrix.Length; i++)
            {
                Console.Write("  ");
                for (int j = 0; j < model.CoefficientMatrix[i].Length; j++)
                {
                    if (j > 0 && model.CoefficientMatrix[i][j] >= 0)
                        Console.Write(" + ");
                    else if (model.CoefficientMatrix[i][j] < 0)
                        Console.Write(" ");
                        
                    Console.Write($"{model.CoefficientMatrix[i][j]:F2}x{j + 1}");
                }
                
                string opStr = model.ConstraintTypes[i] switch
                {
                    ConstraintType.LessThanOrEqual => "<=",
                    ConstraintType.GreaterThanOrEqual => ">=",
                    ConstraintType.Equal => "=",
                    _ => "?"
                };
                
                Console.WriteLine($" {opStr} {model.RHSVector[i]:F2}");
            }
            
            // Display variable constraints
            Console.WriteLine("Where:");
            for (int j = 0; j < model.ObjectiveCoefficients.Length; j++)
            {
                Console.WriteLine($"  x{j + 1} >= 0");
            }
        }
        
        #endregion
    }
}

