using System;
using System.Collections.Generic;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Represents a single iteration of the simplex algorithm
    /// </summary>
    public interface ISimplexIteration
    {
        /// <summary>
        /// Gets the iteration number
        /// </summary>
        int Iteration { get; }
        
        /// <summary>
        /// Gets the current tableau (may be null for revised simplex)
        /// </summary>
        double[,] Tableau { get; }
        
        /// <summary>
        /// Gets the index of the entering variable (-1 if none)
        /// </summary>
        int EnteringVariable { get; }
        
        /// <summary>
        /// Gets the index of the leaving variable (-1 if none)
        /// </summary>
        int LeavingVariable { get; }
        
        /// <summary>
        /// Gets the current objective value
        /// </summary>
        double ObjectiveValue { get; }
        
        /// <summary>
        /// Gets the current solution vector
        /// </summary>
        double[] Solution { get; }
        
        /// <summary>
        /// Gets the current basis variables
        /// </summary>
        int[] Basis { get; }
        
        /// <summary>
        /// Gets the names of the basis variables
        /// </summary>
        IList<string> BasisVariables { get; }
        
        /// <summary>
        /// Gets the names of the non-basis variables
        /// </summary>
        string[] NonBasisVariables { get; }
        
        /// <summary>
        /// Gets the names of all variables
        /// </summary>
        string[] VariableNames { get; }
        
        /// <summary>
        /// Gets the inverse of the basis matrix
        /// </summary>
        double[,] BasisInverse { get; }
        
        /// <summary>
        /// Gets whether this iteration is optimal
        /// </summary>
        bool IsOptimal { get; }
        
        /// <summary>
        /// Gets whether the problem is unbounded at this iteration
        /// </summary>
        bool IsUnbounded { get; }
        
        /// <summary>
        /// Gets whether the problem is infeasible at this iteration
        /// </summary>
        bool IsInfeasible { get; }
        
        /// <summary>
        /// Gets the reduced costs for all variables
        /// </summary>
        double[] ReducedCosts { get; }
        
        /// <summary>
        /// Gets the simplex multipliers (shadow prices)
        /// </summary>
        double[] SimplexMultipliers { get; }
        
        /// <summary>
        /// Gets a description of the iteration
        /// </summary>
        string Description { get; }
        
        /// <summary>
        /// Gets the pivot element used in this iteration
        /// </summary>
        double PivotElement { get; }
    }
}
