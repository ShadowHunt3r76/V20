using System;
using System.Collections.Generic;

namespace LinearProgramming.Algorithms.SensitivityAnalysis
{
    /// <summary>
    /// Defines the standard interface for sensitivity analysis across all solvers
    /// </summary>
    public interface ISensitivityAnalysis
    {
        /// <summary>
        /// Performs the complete sensitivity analysis
        /// </summary>
        void PerformAnalysis();

        /// <summary>
        /// Gets the ranges for a given constraint's RHS
        /// </summary>
        /// <param name="constraintIndex">Index of the constraint</param>
        /// <returns>Tuple of (lowerBound, upperBound) for the RHS</returns>
        (double lowerBound, double upperBound) GetRHSRange(int constraintIndex);

        /// <summary>
        /// Gets the ranges for a given variable's coefficient in the objective function
        /// </summary>
        /// <param name="variableIndex">Index of the variable</param>
        /// <returns>Tuple of (lowerBound, upperBound) for the coefficient</returns>
        (double lowerBound, double upperBound) GetObjectiveCoefficientRange(int variableIndex);

        /// <summary>
        /// Gets the reduced cost for a given variable
        /// </summary>
        /// <param name="variableIndex">Index of the variable</param>
        /// <returns>Reduced cost value</returns>
        double GetReducedCost(int variableIndex);

        /// <summary>
        /// Gets the shadow price for a given constraint
        /// </summary>
        /// <param name="constraintIndex">Index of the constraint</param>
        /// <returns>Shadow price value</returns>
        double GetShadowPrice(int constraintIndex);

        /// <summary>
        /// Gets the current solution values for all variables
        /// </summary>
        /// <returns>Dictionary mapping variable names to their values</returns>
        Dictionary<string, double> GetSolutionValues();

        /// <summary>
        /// Generates a visualization of the sensitivity analysis results
        /// </summary>
        /// <returns>String containing the visualization (e.g., formatted text, chart, etc.)</returns>
        string GenerateVisualization();
    }
}
