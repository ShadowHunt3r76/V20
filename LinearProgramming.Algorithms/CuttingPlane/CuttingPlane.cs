using LinearProgramming.Algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CuttingPlaneAlgorithm
{
    public class CuttingPlane
    {


        public LinearProgramSolution CuttingPlaneSolve(LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel model)
        {
            //Solve LP relaxation using primal simplex
            var simplex = new LinearProgramming.Algorithms.PrimalSimplexSolver();
            var solution = simplex.Solve(model);
            double[,] tableau = solution.OptimalTable;

            int numOriginalVars = model.ObjectiveCoefficients.Length;

            //Save initial tableau (LP relaxation)
            SaveIteration(tableau);

            bool integerSolutionFound = false;
            int maxIter = 250;
            int iter = 0;

            //Cutting plane loop
            while (!integerSolutionFound && iter < maxIter)
            {
                //Find fractional row among decision variables
                int cutRow = SelectPivotRowCuttingP(tableau, numOriginalVars);
                if (cutRow == -1)
                {
                    integerSolutionFound = true;
                    break;
                }

                //Apply cut
                tableau = ApplyCuttingP(tableau, cutRow);
                SaveIteration(tableau);

                //Restore feasibility with Dual Simplex
                bool feasible = false;
                int dualIter = 0;
                while (!feasible && dualIter < 50)
                {
                    feasible = OptimalityCheckDualS(tableau);
                    SaveIteration(tableau);
                    dualIter++;
                }

                iter++;
            }

            //Print final results
            PrintAnswersToTxt("cuttingplane_output.txt");

            return new LinearProgramSolution
            {
                Status = integerSolutionFound ? "Optimal Integer Solution (Cutting Plane)" : "Max iterations reached",
                SolutionVector = ExtractSolution(tableau, numOriginalVars),
                ObjectiveValue = tableau[0, tableau.GetLength(1) - 1],
                OptimalTable = tableau,
                TableHistory = CuttingPlaneHistory
            };
        }

        private double[] ExtractSolution(double[,] tableau, int numOriginalVars)
        {
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
                    solution[j] = tableau[basicRow, rhsIndex];
                else
                    solution[j] = 0.0;
            }

            return solution;
        }

        //cutting plane pivot row selecection
        private int SelectPivotRowCuttingP(double[,] tableau, int numOriginalVariables)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            int rhsIndex = cols - 1;

            int bestRow = -1;
            double bestFraction = 0.0;

            // loop over all constraint rows
            for (int i = 1; i < rows; i++)
            {
                double rhs = tableau[i, rhsIndex];
                double fraction = rhs - Math.Floor(rhs);

                if (fraction > 1e-6) // not an integer
                {
                    // find the basic variable in this row
                    int basicCol = -1;
                    for (int j = 0; j < cols - 1; j++)
                    {
                        bool isUnitColumn = true;
                        if (tableau[i, j] == 1)
                        {
                            // check other rows for 0
                            for (int k = 0; k < rows; k++)
                            {
                                if (k != i && Math.Abs(tableau[k, j]) > 1e-6)
                                {
                                    isUnitColumn = false;
                                    break;
                                }
                            }
                            if (isUnitColumn)
                            {
                                basicCol = j;
                                break;
                            }
                        }
                    }

                    // only allow if the basic variable is an original x-variable
                    if (basicCol != -1 && basicCol < numOriginalVariables)
                    {
                        if (fraction > bestFraction)
                        {
                            bestFraction = fraction;
                            bestRow = i;
                        }
                    }
                }
            }

            return bestRow; // -1 means "no fractional x-variable found"
        }


        //cutting plane cut method

        private double[,] ApplyCuttingP(double[,] tableau, int cutRow)
        {
            // Helper: fractional part function
            double FractionalPart(double x)
            {
                double f = x - Math.Floor(x);   // always between 0 and 1
                if (Math.Abs(f) < 1e-6) return 0;  // treat tiny value as 0
                return f;
            }

            int oldRows = tableau.GetLength(0);
            int oldCols = tableau.GetLength(1);
            int rhsIndex = oldCols - 1;

            // New tableau: +1 row for the cut, +1 column for the slack variable
            double[,] newTableau = new double[oldRows + 1, oldCols + 1];

            // Copy the old tableau (coefficients)
            for (int i = 0; i < oldRows; i++)
            {
                for (int j = 0; j < oldCols; j++)
                {
                    newTableau[i, j] = tableau[i, j];
                }
            }

            // Copy RHS values into new last column
            for (int i = 0; i < oldRows; i++)
            {
                newTableau[i, oldCols] = tableau[i, rhsIndex];
            }

            // Build the cut row
            for (int j = 0; j < rhsIndex; j++)
            {
                double coeff = tableau[cutRow, j];
                newTableau[oldRows, j] = -FractionalPart(coeff); // negative fractional part
            }

            // Slack variable coefficient = 1
            newTableau[oldRows, oldCols - 1] = 1;

            // RHS = - fractional part of original RHS
            double rhs = tableau[cutRow, rhsIndex];
            newTableau[oldRows, oldCols] = -FractionalPart(rhs);

            return newTableau;
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

        private bool OptimalityCheckDualS(double[,] tableau)
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

            // If we reach here → all RHS are ≥ 0
            Console.WriteLine("Dual Simplex: Feasible solution found.");
            // Here you would normally call OptimalityCheckPrimal(tableau) and maybe PrintAnswersToTxt()
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



        //print answers to txt file

        // Store history of cutting plane iterations
        private List<double[,]> CuttingPlaneHistory = new List<double[,]>();

        private void SaveIteration(double[,] tableau)
        {
            // clone the tableau so later changes don't overwrite
            CuttingPlaneHistory.Add((double[,])tableau.Clone());
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
