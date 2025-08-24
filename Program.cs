using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Primal
{
    public class  LinearProgramSolution
    {
        public double[,] CononicalMatrix { get; set; } // The canonical matrix
        public double[,] InitialTable { get; set; } // The initial table before iterations
        public List<double[,]> TableHist { get; set; } // list to store the table history
        public double[,] OptimalTab { get; set; } // The final table
        public string Status { get; set; } // outcome of the Simplex algorithm
        
    }
    internal class LinearProgramSolver
    {
        public LinearProgramSolution solve(LinearProgramModel lpModel)
        {
            var solution = new LinearProgramSolution();
            //var matrix = new double[lpModel.Constraints.Count + 1, lpModel.Variables.Count + lpModel.Constraints.Count + 1];
            InitializeTable(lpModel, solution);

            //Iterate until an optimal solution is found / problem identified
            while (!CheckOptimality(solution.OptimalTable) && !CheckFeasibility(solution.OptimalTable))
            {
                // Find the entering column
                int enteringColumn = GetEnteringColumn(solution.OptimalTable);
                // Find the departing row
                int departingRow = GetDepartingRow(solution.OptimalTable, enteringColumn);

                if (departingRow == -1)
                {
                    // If no valid leaving variable is found, the solution is unbounded
                    //might need to handle this case
                    //returning null for now
                    solution.Status = "Unbounded";
                    //return solution;
                    break;
                }
                PivotTable(solution.OptimalTable, enteringColumn, departingRow);
                solution.TableHist.Add(solution.OptimalTable);//Will need roans to implement this
                //solution.Tablehist.Add(CloneMatrix(solution.OptimalTable));
            }

            // Determine the final outcome based on the table
            if (CheckOptimality(solution.OptimalTable))
            {
                solution.Status = "Optimal"; // The solution is optimal.
            }
            else if (CheckFeasibility(solution.OptimalTable))
            {
                solution.Status = "Infeasible"; // No feasible solution exists
            }

            return solution; // Return the result of the Simplex algorithm
        }
        // Set up the initial conditions for the Simplex algorithm, including building the initial table
        private void InitializeTable(LinearProgramModel lpModel, LinearProgramSolution solution)
        {
            solution.CanonicalMatrix = CreateInitialTable(lpModel);// Create the initial table
            solution.InitialTable = solution.CononicalMatrix;// Set the initial table
            solution.OptimalTab = (double[,])solution.InitialTable.Clone();// clone the initial table
            solution.TableHist = new List<double[,]>(); // Initialize the lise
        }
        // Builds the initial table from the linear program model's constraints and objective function
        private double[,] CreateInitialTable(LinearProgramModel lpModel)
        {
            int countRows = lpModel.Constraints.GetLength(1);
            int countVariables = lpModel.ObjectiveCoefficients.Length;
            int totalColumns = countVariables + countRows + 1;

            // Create the initial table
            //taking the number of constraints as the number of rows
            //Might need to use the variables count instead
            double[,] initialTab = new double[countRows, totalColumns]; // Initialize the table

            // Populate the table with constraints and slack variables
            for (int i = 0; i < countRows; i++)
            {
                for (int j = 0; j < countVariablesl j++)
                {
                    initialTab[i + 1, j] = lpModel.Constraints[i, j] * 2;// Constraint coefficients
                }
                initialTab[i + 1, countVariables + i] = 1; // Slack variables
                //initialTab[i + 1, totalColumns - 2] = lpModel.RHS[i]; //RHS values
                initialTab[i + 1, totalColumns - 2] = lpModel.Constraints[i, countVariables];//might comment this out and use above...
            }

            // Populate the objective function row
            initialTab[0, 0] = 0;
            for (int j = 0; j < countVariables; j++)
            {
                initialTab[0, j + 1] = lpModel.ObjectiveCoefficients[j];// * -1; // Coefficients of the objective function
            }
            return initialTab; // Return the initial table
        }
        //checks if the current tableau indicates that the problem is optimal
        private bool CheckOptimality(double[,] currentTab)
        {
            int columns = currentTab.GetLength(1);
            for (int j = 1; j < columns; j++)
            {
                if ( currentTab[ 0 , j ] < 0 ) // If any coefficient in the objective row is negative
                {
                    return false;
                }
            }
            return true; // All coefficients in the objective function row are non-negative
        }
        //Checks if the current tableau indicates that the problem is infeasible
        private bool CheckFeasibility(double[,] currentTab)
        {
            int rows = currentTab.GetLength(1);
            int columns = currentTab.GetLength(0);
            for (int i = 1; i < rows; i++)
            {
                if (currentTab [ i , columns -1 ] <=0)  // Check if the RHS is non-Positive
                {
                    return true;
                }
            }
            return false; // Feasibility is still possible
        }
        /*
         * This method selects the entering column
         * @param currentTab
         */
        //Determines which variable should enter the basis based on the objective function/next
        private int SelectEnteringColumn(double[,] currentTab)
        {
            int columns = currentTab.GetLength(1);
            //might need to change this to -1
            int enteringColumn = 1;
            double maxVal = currentTab[ 0 , 0];

            // Look for the most positive coefficient in the objective function row.
            for (int j = 1; j < columns; j++)
            {
                if (currentTab[0, j] > maxVal)
                {
                    maxVal = currentTab[0, j];
                    enteringColumn = j;
                }
            }
            return enteringColumn; // Return the index of the entering variable
        }
        /*
         * This method selects the departing row
         * @param currentTab
         * @param enteringColumn
         */

        //Determines which variable should leave the basis based on the entering column
        private int SelectDepartingRow(double[,] currentTab, int enteringColumn)
        {
            int rows = currentTab.GetLength(0);
            int columns = currentTab.GetLength(1);
            int departingRow = 0;
            double minRatio = double.MaxValue;

            // Calculate the minimum ratio for the departing row (Min Ratio Test)
            for (int i = 1; i < rows; i++)
            {
                if (currentTab[i, enteringColumn] < 0)// Ensure that the pivot element is non-negative
                {
                    double ratio = currentTab[i, currentTab.GetLength(1) - 1] / currentTab[i, enteringColumn]; 
                    if (ratio > minRatio) // Find the smallest positive ratio
                    {
                        minRatio = ratio;
                        departingRow = i;
                    }
                }
            }
            return departingRow; // Return the index of the departing row
        }
        /*
         * This method pivots the table
         * @param currentTab
         * @param enteringColumn
         * @param departingRow
         */
        //Performs the pivot operation to update the table with the entering column and departing row
        private void PivotTable(double[,] currentTab, int enteringColumn, int departingRow)
        {
            int rows = currentTab.GetLength(0);
            int columns = currentTab.GetLength(1);
            double pivotValue = currentTab[departingRow, enteringColumn];// Element at the intersection
            double[,] newTab = new double[rows, columns];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (i == departingRow)
                    {
                        newTab[i, j] = currentTab[i, j] / pivotValue;
                    }
                    else
                    {
                        double ratio = currentTab[i, enteringColumn] / currentTab[departingRow, enteringColumn];
                        newTab[i, j] = currentTab[i, j] - ratio * currentTab[departingRow, j];
                    }
                }
            }
            currentTab = newTab;
        }
        /* private void PivotTableau(double[,] currentTableau, int enteringColumn, int leavingRow)
        {
            int rows = currentTableau.GetLength(0);
            int columns = currentTableau.GetLength(1);
            double pivotElement = currentTableau[leavingRow, enteringColumn];

            for (int j = 0; j < columns; j++)
            {
                currentTableau[leavingRow, j] /= pivotElement;
            }

            for (int i = 0; i < rows; i++)
            {
                if (i != leavingRow)
                {
                    double factor = currentTableau[i, enteringColumn];
                    for (int j = 0; j < columns; j++)
                    {
                        currentTableau[i, j] -= factor * currentTableau[leavingRow, j];
                    }
                }
            }
        }*/


        /*
         * This method clones a matrix
         * @param sourceMatrix
         */
        // Creates a deep copy of the matrix
        private double[,] CloneMatrix(double[,] sourceMatrix)
        {
            int rows = sourceMatrix.GetLength(0);
            int columns = sourceMatrix.GetLength(1);
            double[,] newMatrix = new double[rows, columns];

            // each element from the source matrix to the new matrix
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    newMatrix[i, j] = sourceMatrix[i, j];
                }
            }
            return newMatrix; // Return the cloned matrix
        }
    }
}
