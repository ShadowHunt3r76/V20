using System;
using System.Collections.Generic;
using LinearProgramming.Parsing;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Represents the solution to a linear programming problem
    /// </summary>
    public class LinearProgramSolution : ILinearProgramSolution
    {
        /// <summary>
        /// Gets or sets the status of the solution (e.g., "Optimal", "Unbounded", "Infeasible")
        /// </summary>
        public string Status { get; set; }

        /// <summary>
        /// Returns a detailed string representation of the solution
        /// </summary>
        public override string ToString()
        {
            var sb = new System.Text.StringBuilder();
            
            sb.AppendLine("\n=== Linear Programming Solution ===");
            sb.AppendLine($"Status: {Status}");
            sb.AppendLine($"Optimal Value: {ObjectiveValue:F4}\n");
            
            if (VariableValues != null && VariableValues.Count > 0)
            {
                sb.AppendLine("Variable Values:");
                foreach (var kvp in VariableValues)
                {
                    sb.AppendLine($"  {kvp.Key} = {kvp.Value:F4}");
                }
                sb.AppendLine();
            }
            
            if (BasisVariables != null && BasisVariables.Length > 0)
            {
                sb.AppendLine("Basis Variables:");
                for (int i = 0; i < BasisVariables.Length; i++)
                {
                    sb.AppendLine($"  x{i + 1} = {BasisVariables[i]}");
                }
                sb.AppendLine();
            }
            
            return sb.ToString();
        }
        
        /// <summary>
        /// Displays the complete solution including all iterations
        /// </summary>
        public void DisplaySolution()
        {
            Console.WriteLine("\n=== Linear Programming Solution ===");
            Console.WriteLine($"Status: {Status}");
            Console.WriteLine($"Optimal Value: {ObjectiveValue:F4}\n");
            
            if (Iterations != null && Iterations.Count > 0)
            {
                Console.WriteLine($"Total Iterations: {Iterations.Count}\n");
                
                // Display each iteration
                foreach (var iteration in Iterations)
                {
                    Console.WriteLine($"=== Iteration {iteration.Iteration} ===");
                    Console.WriteLine($"Status: {iteration.Description}");
                    
                    if (iteration.EnteringVariable >= 0 && VariableNames != null && 
                        iteration.EnteringVariable < VariableNames.Length)
                    {
                        Console.WriteLine($"Entering: {VariableNames[iteration.EnteringVariable]} (x{iteration.EnteringVariable + 1})");
                    }
                    
                    if (iteration.LeavingVariable >= 0 && BasisVariables != null && 
                        iteration.LeavingVariable < BasisVariables.Length)
                    {
                        Console.WriteLine($"Leaving: {BasisVariables[iteration.LeavingVariable]} (row {iteration.LeavingVariable + 1})");
                    }
                    
                    Console.WriteLine($"Objective Value: {iteration.ObjectiveValue:F4}\n");
                }
            }
            
            // Display final solution
            Console.WriteLine("\n=== Final Solution ===");
            Console.WriteLine(ToString());
        }

        /// <summary>
        /// Gets or sets the optimal value of the objective function
        /// </summary>
        public double ObjectiveValue { get; set; }

        /// <summary>
        /// Gets or sets the values of the decision variables in the optimal solution
        /// </summary>
        public Dictionary<string, double> VariableValues { get; set; } = new Dictionary<string, double>();

        /// <summary>
        /// Gets or sets the basis indices in the final solution
        /// </summary>
        public int[] Basis { get; set; }

        /// <summary>
        /// Gets or sets the reduced costs in the final solution
        /// </summary>
        public double[] ReducedCosts { get; set; }

        /// <summary>
        /// Gets or sets the simplex multipliers (dual variables)
        /// </summary>
        public double[] SimplexMultipliers { get; set; }

        /// <summary>
        /// Gets or sets the list of simplex iterations
        /// </summary>
        public List<ISimplexIteration> Iterations { get; set; } = new List<ISimplexIteration>();

        /// <summary>
        /// Gets or sets the list of detailed iteration information for revised simplex
        /// </summary>
        public List<RevisedSimplexIteration> IterationDetails { get; set; } = new List<RevisedSimplexIteration>();

        /// <summary>
        /// Gets the number of iterations performed
        /// </summary>
        public int IterationCount => Iterations?.Count ?? 0;

        /// <summary>
        /// Gets or sets the solution vector containing the values of all variables
        /// </summary>
        public double[] SolutionVector { get; set; }

        /// <summary>
        /// Gets or sets the initial simplex tableau
        /// </summary>
        public double[,] InitialTable { get; set; }

        /// <summary>
        /// Gets or sets the final simplex tableau
        /// </summary>
        public double[,] OptimalTable { get; set; }

        /// <summary>
        /// Gets or sets the history of tableaus during the simplex algorithm
        /// </summary>
        public List<double[,]> TableHistory { get; set; } = new List<double[,]>();
        
        /// <summary>
        /// Gets or sets the names of the basis variables in the final tableau
        /// </summary>
        public string[] BasisVariables { get; set; }
        
        /// <summary>
        /// Gets or sets the names of all variables (including slacks and artificials)
        /// </summary>
        public string[] VariableNames { get; set; }

        /// <summary>
        /// Gets or sets the canonical matrix (A) of the linear program
        /// </summary>
        public double[,] CanonicalMatrix { get; set; }
        
        /// <summary>
        /// Gets or sets the names of non-basis variables
        /// </summary>
        public string[] NonBasisVariables { get; set; }

        /// <summary>
        /// Gets a value indicating whether the solution is optimal
        /// </summary>
        public bool IsOptimal => Status == "Optimal";

        /// <summary>
        /// Gets a value indicating whether the problem is unbounded
        /// </summary>
        public bool IsUnbounded => Status == "Unbounded";

        /// <summary>
        /// Gets a value indicating whether the problem is infeasible
        /// </summary>
        public bool IsInfeasible => Status == "Infeasible";

        /// <summary>
        /// Gets or sets a value indicating whether the solution is degenerate
        /// </summary>
        public bool IsDegenerate { get; set; }


        /// <summary>
        /// Gets or sets a value indicating whether to use the product form of the inverse
        /// </summary>
        public bool UsingProductForm { get; set; } = true;
        
        /// <summary>
        /// Gets or sets the basis inverse matrix
        /// </summary>
        public double[,] BasisInverse { get; set; }
    }
}
