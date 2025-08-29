using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Parsing;

namespace LinearProgramming.Algorithms
{
    // Use the outer CanonicalLinearProgrammingModel from LinearProgramming.Parsing
    using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;

    /// <summary>
    /// Implements the Revised Primal Simplex algorithm for solving linear programs.
    /// This implementation includes proper matrix operations, basis management, and Phase I/II method.
    /// </summary>
    public class RevisedPrimalSimplexSolver
    {
        private const double EPSILON = 1e-10;
        private const double BIG_M = 1e6;

        public LinearProgramSolution Solve(CanonicalLinearProgrammingModel model)
        {
            int m = model.CoefficientMatrix.Length; // Number of constraints
            int n = model.CoefficientMatrix[0].Length; // Number of variables
            double[][] A = model.CoefficientMatrix;
            double[] b = model.RHSVector;
            double[] c = model.ObjectiveCoefficients;
            var variableTypes = model.VariableTypes;
            
            var solution = new LinearProgramSolution();
            
            // Check if any RHS values are negative - convert to positive if needed
            for (int i = 0; i < m; i++)
            {
                if (b[i] < 0)
                {
                    b[i] = -b[i];
                    for (int j = 0; j < n; j++)
                    {
                        A[i][j] = -A[i][j];
                    }
                }
            }
            
            // Find initial feasible basis - first try to identify natural basis
            var initialBasisResult = FindInitialBasis(A, b, variableTypes, m, n);
            List<int> basis = initialBasisResult.Basis;
            bool needsPhaseI = initialBasisResult.NeedsPhaseI;
            double[][] workingA = initialBasisResult.WorkingA;
            double[] workingC = initialBasisResult.WorkingC;
            double[] workingB = initialBasisResult.WorkingB;
            
            // Phase I if artificial variables are needed
            if (needsPhaseI)
            {
                var phaseIResult = SolvePhaseI(workingA, workingB, basis, m);
                if (phaseIResult.Status != "Optimal")
                {
                    solution.Status = phaseIResult.Status;
                    return solution;
                }
                
                // Check if artificial variables are zero
                bool hasNonZeroArtificial = false;
                for (int i = 0; i < basis.Count; i++)
                {
                    if (basis[i] >= n && phaseIResult.BasicSolution[i] > EPSILON)
                    {
                        hasNonZeroArtificial = true;
                        break;
                    }
                }
                
                if (hasNonZeroArtificial)
                {
                    solution.Status = "Infeasible";
                    return solution;
                }
                
                // Remove artificial variables and continue to Phase II
                basis = phaseIResult.Basis;
                workingB = phaseIResult.BasicSolution;
            }
            
            // Phase II - solve original problem
            return SolvePhaseII(A, b, c, basis, m, n);
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
        
        /// <summary>
        /// Dot product of two vectors
        /// </summary>
        private double DotProduct(double[] a, double[] b)
        {
            double result = 0;
            for (int i = 0; i < a.Length; i++)
            {
                result += a[i] * b[i];
            }
            return result;
        }
        
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
    }
}
