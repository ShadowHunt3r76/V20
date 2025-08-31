using System;
using System.Collections.Generic;

namespace LinearProgramming.Parsing
{
    /// <summary>
    /// Represents the solution to a linear programming problem
    /// </summary>
    public interface ILinearProgramSolution
    {
        /// <summary>
        /// Gets or sets the status of the solution (e.g., "Optimal", "Unbounded", "Infeasible")
        /// </summary>
        string Status { get; set; }

        /// <summary>
        /// Gets or sets the optimal value of the objective function
        /// </summary>
        double ObjectiveValue { get; set; }

        /// <summary>
        /// Gets or sets the values of the decision variables in the optimal solution
        /// </summary>
        Dictionary<string, double> VariableValues { get; set; }

        /// <summary>
        /// Gets or sets the solution vector containing the values of all variables
        /// </summary>
        double[] SolutionVector { get; set; }

        /// <summary>
        /// Gets or sets the initial tableau
        /// </summary>
        double[,] InitialTable { get; set; }

        /// <summary>
        /// Gets or sets the optimal tableau
        /// </summary>
        double[,] OptimalTable { get; set; }

        /// <summary>
        /// Gets or sets the history of tableaus during optimization
        /// </summary>
        List<double[,]> TableHistory { get; set; }

        /// <summary>
        /// Gets or sets the names of the basis variables in the final tableau
        /// </summary>
        string[] BasisVariables { get; set; }

        /// <summary>
        /// Gets or sets the canonical matrix (A) of the linear program
        /// </summary>
        double[,] CanonicalMatrix { get; set; }

        /// <summary>
        /// Gets or sets the names of all variables (including slacks and artificials)
        /// </summary>
        string[] VariableNames { get; set; }
    }
}
