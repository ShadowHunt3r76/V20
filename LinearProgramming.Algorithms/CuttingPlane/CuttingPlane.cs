using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CuttingPlaneAlgorithm
{
    internal class CuttingPlane
    {
        

        public CuttingPlane() { }

        //cutting plane pivot row selecection
        public void SelectPivotRowCuttingP(double[,] tableau, int[] basicVariables, int[] nonBasicVariables)
        {
            // Implement the logic to select the pivot row based on the cutting plane method
            // This typically involves finding the most negative coefficient in the objective function row
            // and determining which variable to enter the basis.
            // For now, this is a placeholder for your implementation.
        }

        //cutting plane cut method

        public void ApplyCuttingP(double[,] tableau, int[] basicVariables, int[] nonBasicVariables)
        {
            // Implement the logic to apply the cutting plane method
            // This typically involves adding a new constraint to the tableau based on the current solution
            // and resolving the tableau.
            // Add new constraint to optimal table.
        }

        //Dual Simplex programming

        // pivot row selection DualSimplex

        public void SelectPivotRowDualS(double[,] tableau, int[] basicVariables, int[] nonBasicVariables)
        {
            // Implement the logic to select the pivot row for the dual simplex method
            // This typically involves finding the most negative coefficient in the right-hand side
            // and determining which variable to enter the basis.
            // For now, this is a placeholder for your implementation.

            //call SelectPivotColumnDualS(tableau, basicVariables, nonBasicVariables);
        }

        //pivot column selection DualSimplex

        public void SelectPivotColumnDualS(double[,] tableau, int[] basicVariables, int[] nonBasicVariables)
        {
            // Implement the logic to select the pivot column for the dual simplex method
            // This typically involves finding the most negative coefficient in the objective function row
            // and determining which variable to leave the basis.
            // For now, this is a placeholder for your implementation.

            //call PivotPrimalSimplex method;
        }

        //PivotPrimaSimplex method call here & call OptimalityCheckDualS method

        // optimality check Dual Simplex

        public void OptimalityCheckDualS(double[,] tableau, int[] basicVariables, int[] nonBasicVariables)
        {
            // Implement the logic to check the optimality of the current solution in the dual simplex method
            // This typically involves checking if all coefficients in the objective function row are non-negative.
            // For now, this is a placeholder for your implementation.

            //If optimal for Dual, call OptimalityCheckPrimal; If optimal - PrintAnswersToTxt method call.
            //If not optimal for Dual, call SelectPivotRowDualS method;
        }

        //print answers to txt file

        public void PrintAnswersToTxt()
        {
            
        }




    }
}
