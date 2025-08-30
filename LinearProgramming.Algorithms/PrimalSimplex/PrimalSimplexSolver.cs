using System;
using System.Collections.Generic;
using LinearProgramming.Parsing;

namespace LinearProgramming.Algorithms
{
    // Use the outer CanonicalLinearProgrammingModel from LinearProgramming.Parsing
    using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;

    public class LinearProgramSolution
    {
        public double[,] CanonicalMatrix { get; set; }
        public double[,] InitialTable { get; set; }
        public List<double[,]> TableHistory { get; set; } = new List<double[,]>();
        public double[,] OptimalTable { get; set; }
        public string Status { get; set; }
        public double[] SolutionVector { get; set; }
        public double ObjectiveValue { get; set; }
    }

    public class PrimalSimplexSolver
    {
        public LinearProgramSolution Solve(CanonicalLinearProgrammingModel model)
        {
            int m = model.CoefficientMatrix.Length;
            int n = model.CoefficientMatrix[0].Length;
            // Build initial tableau
            double[,] tableau = BuildInitialTableau(model);
            var solution = new LinearProgramSolution
            {
                CanonicalMatrix = To2D(model.CoefficientMatrix),
                InitialTable = (double[,])tableau.Clone(),
                TableHistory = new List<double[,]>()
            };
            int maxIter = 1000;
            int iter = 0;
            while (!IsOptimal(tableau) && iter < maxIter)
            {
                int entering = GetEnteringColumn(tableau);
                int leaving = GetLeavingRow(tableau, entering);
                if (leaving == -1)
                {
                    solution.Status = "Unbounded";
                    solution.OptimalTable = tableau;
                    return solution;
                }
                tableau = Pivot(tableau, entering, leaving);
                solution.TableHistory.Add((double[,])tableau.Clone());
                iter++;
            }
            solution.OptimalTable = tableau;
            solution.Status = IsOptimal(tableau) ? "Optimal" : "Infeasible";
            solution.SolutionVector = ExtractSolution(tableau);
            solution.ObjectiveValue = tableau[0, tableau.GetLength(1) - 1];
            return solution;
        }

        private double[,] BuildInitialTableau(CanonicalLinearProgrammingModel model)
        {
            int m = model.CoefficientMatrix.Length;
            int n = model.CoefficientMatrix[0].Length;
            double[,] tableau = new double[m + 1, n + 1];
            // Objective row
            tableau[0, 0] = 0;
            for (int j = 0; j < n; j++)
                tableau[0, j] = -model.ObjectiveCoefficients[j];
            // Constraints
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                    tableau[i + 1, j] = model.CoefficientMatrix[i][j];
                tableau[i + 1, n] = model.RHSVector[i];
            }
            return tableau;
        }

        private bool IsOptimal(double[,] tableau)
        {
            int cols = tableau.GetLength(1);
            for (int j = 0; j < cols - 1; j++)
                if (tableau[0, j] < 0) return false;
            return true;
        }

        private int GetEnteringColumn(double[,] tableau)
        {
            int cols = tableau.GetLength(1);
            int idx = -1;
            double min = 0;
            for (int j = 0; j < cols - 1; j++)
            {
                if (tableau[0, j] < min)
                {
                    min = tableau[0, j];
                    idx = j;
                }
            }
            return idx;
        }

        private int GetLeavingRow(double[,] tableau, int entering)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            int idx = -1;
            double minRatio = double.MaxValue;
            for (int i = 1; i < rows; i++)
            {
                double coeff = tableau[i, entering];
                if (coeff > 0)
                {
                    double ratio = tableau[i, cols - 1] / coeff;
                    if (ratio < minRatio)
                    {
                        minRatio = ratio;
                        idx = i;
                    }
                }
            }
            return idx;
        }

        private double[,] Pivot(double[,] tableau, int entering, int leaving)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            double[,] newTab = new double[rows, cols];
            double pivot = tableau[leaving, entering];
            for (int j = 0; j < cols; j++)
                newTab[leaving, j] = tableau[leaving, j] / pivot;
            for (int i = 0; i < rows; i++)
            {
                if (i == leaving) continue;
                double factor = tableau[i, entering];
                for (int j = 0; j < cols; j++)
                    newTab[i, j] = tableau[i, j] - factor * newTab[leaving, j];
            }
            return newTab;
        }

        //made public, so code is reusable; implemented in cutting plane
        public double[] ExtractSolution(double[,] tableau)
        {
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            double[] solution = new double[cols - 1];
            for (int j = 0; j < cols - 1; j++)
            {
                solution[j] = 0;
                for (int i = 1; i < rows; i++)
                {
                    bool isBasic = true;
                    for (int k = 1; k < rows; k++)
                        if (k != i && tableau[k, j] != 0) isBasic = false;
                    if (isBasic && tableau[i, j] == 1)
                    {
                        solution[j] = tableau[i, cols - 1];
                        break;
                    }
                }
            }
            return solution;
        }

        private double[,] To2D(double[][] arr)
        {
            int rows = arr.Length;
            int cols = arr[0].Length;
            double[,] result = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    result[i, j] = arr[i][j];
            return result;
        }
    }
}
