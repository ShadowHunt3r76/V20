using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.Utils;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Represents the solution to a linear programming problem
    /// </summary>
    public class LinearProgramSolution : ILinearProgramSolution
    {
        #region ILinearProgramSolution Implementation
        
        /// <summary>
        /// Gets or sets the status of the solution (e.g., "Optimal", "Unbounded", "Infeasible")
        /// </summary>
        public string Status { get; set; } = string.Empty;

        /// <summary>
        /// Gets or sets the optimal value of the objective function
        /// </summary>
        public double ObjectiveValue { get; set; }

        /// <summary>
        /// Gets or sets the values of the decision variables in the optimal solution
        /// </summary>
        public Dictionary<string, double> VariableValues { get; set; } = new();

        /// <summary>
        /// Gets or sets the solution vector containing the values of all variables
        /// </summary>
        public double[] SolutionVector { get; set; } = Array.Empty<double>();

        /// <summary>
        /// Gets or sets the initial tableau
        /// </summary>
        public double[,] InitialTable { get; set; } = new double[0, 0];

        /// <summary>
        /// Gets or sets the optimal tableau
        /// </summary>
        public double[,] OptimalTable { get; set; } = new double[0, 0];

        /// <summary>
        /// Gets or sets the history of tableaus during optimization
        /// </summary>
        public List<double[,]> TableHistory { get; set; } = new();

        /// <summary>
        /// Gets or sets the names of the basis variables in the final tableau
        /// </summary>
        public string[] BasisVariables { get; set; } = Array.Empty<string>();

        /// <summary>
        /// Gets or sets the indices of the basis variables in the final tableau
        /// </summary>
        public int[] BasisIndices { get; set; } = Array.Empty<int>();

        /// <summary>
        /// Gets or sets the canonical matrix (A) of the linear program
        /// </summary>
        public double[,] CanonicalMatrix { get; set; } = new double[0, 0];

        /// <summary>
        /// Gets or sets the names of all variables (including slacks and artificials)
        /// </summary>
        public string[] VariableNames { get; set; } = Array.Empty<string>();
        
        #endregion
        
        #region Additional Properties
        
        /// <summary>
        /// Gets or sets the objective value of the LP relaxation (used in integer programming)
        /// </summary>
        public double? LPObjectiveValue { get; set; }
        
        /// <summary>
        /// Gets or sets the basis indices in the final solution
        /// </summary>
        public int[] Basis { get; set; }

        /// <summary>
        /// Gets or sets the reduced costs in the final solution
        /// </summary>
        public double[] ReducedCosts { get; set; } = Array.Empty<double>();

        /// <summary>
        /// Gets or sets the simplex multipliers (dual variables)
        /// </summary>
        public double[] SimplexMultipliers { get; set; } = Array.Empty<double>();

        /// <summary>
        /// Gets or sets the list of simplex iterations
        /// </summary>
        public List<ISimplexIteration> Iterations { get; set; } = new List<ISimplexIteration>();
        
        /// <summary>
        /// Gets or sets the list of detailed iteration information for revised simplex
        /// </summary>
        public List<RevisedSimplexIteration> IterationDetails { get; set; } = new List<RevisedSimplexIteration>();

        /// <summary>
        /// Gets or sets the names of non-basis variables
        /// </summary>
        public string[] NonBasisVariables { get; set; } = Array.Empty<string>();
        
        /// <summary>
        /// Gets the number of iterations performed
        /// </summary>
        public int IterationCount => Iterations?.Count ?? 0;
        
        /// <summary>
        /// Gets a value indicating whether the solution is optimal
        /// </summary>
        public bool IsOptimal => Status == "Optimal";
        
        #endregion

        /// <summary>
        /// Returns a detailed string representation of the solution
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();
            
            // Add status and objective value
            sb.AppendLine(OutputFormatter.FormatKeyValue("Status", Status));
            sb.AppendLine(OutputFormatter.FormatKeyValue("Optimal Value", ObjectiveValue.ToString("F4")));
            
            // Add variable values
            if (VariableValues != null && VariableValues.Count > 0)
            {
                sb.AppendLine("\nVariable Values:");
                var varTable = VariableValues.Select(kvp => new object[] { kvp.Key, kvp.Value });
                sb.AppendLine(OutputFormatter.CreateTable(
                    new[] { "Variable", "Value" },
                    varTable
                ));
            }
            
            // Add basis variables
            if (BasisVariables != null && BasisVariables.Length > 0)
            {
                sb.AppendLine("\nBasis Variables:");
                var basisTable = BasisVariables.Select((b, i) => new object[] { $"x{i + 1}", b });
                sb.AppendLine(OutputFormatter.CreateTable(
                    new[] { "Position", "Variable" },
                    basisTable
                ));
            }
            
            // Add reduced costs if available
            if (ReducedCosts != null && ReducedCosts.Length > 0 && VariableNames != null)
            {
                sb.AppendLine("\nReduced Costs:");
                var reducedCostsTable = VariableNames
                    .Select((name, i) => new { Name = name, Index = i })
                    .Where(x => x.Index < ReducedCosts.Length)
                    .Select(x => new object[] { x.Name, ReducedCosts[x.Index] });
                    
                sb.AppendLine(OutputFormatter.CreateTable(
                    new[] { "Variable", "Reduced Cost" },
                    reducedCostsTable
                ));
            }
            
            // Add simplex multipliers (dual variables) if available
            if (SimplexMultipliers != null && SimplexMultipliers.Length > 0)
            {
                sb.AppendLine("\nSimplex Multipliers (Dual Variables):");
                var dualsTable = SimplexMultipliers
                    .Select((value, i) => new object[] { $"Constraint {i + 1}", value });
                    
                sb.AppendLine(OutputFormatter.CreateTable(
                    new[] { "Constraint", "Dual Value" },
                    dualsTable
                ));
            }
            
            return OutputFormatter.CreateBox(sb.ToString(), "LINEAR PROGRAMMING SOLUTION");
        }
        
        /// <summary>
        /// Displays the complete solution including all iterations
        /// </summary>
        public void DisplaySolution()
        {
            // Display solution summary
            Console.WriteLine(ToString());
            
            // Display iteration details if available
            if (Iterations != null && Iterations.Count > 0)
            {
                var iterationSummary = new StringBuilder();
                iterationSummary.AppendLine(OutputFormatter.FormatKeyValue("Total Iterations", Iterations.Count));
                
                // Add iteration summary table
                var iterationRows = Iterations.Select(iter => new object[] {
                    iter.Iteration,
                    iter.Description,
                    iter.EnteringVariable >= 0 && VariableNames != null && iter.EnteringVariable < VariableNames.Length 
                        ? VariableNames[iter.EnteringVariable] 
                        : "-",
                    iter.LeavingVariable >= 0 && BasisVariables != null && iter.LeavingVariable < BasisVariables.Length
                        ? BasisVariables[iter.LeavingVariable]
                        : "-",
                    iter.ObjectiveValue.ToString("F4")
                });
                
                iterationSummary.AppendLine(OutputFormatter.CreateTable(
                    new[] { "Iteration", "Status", "Entering", "Leaving", "Objective" },
                    iterationRows
                ));
                
                Console.WriteLine(OutputFormatter.CreateBox(iterationSummary.ToString(), "ITERATION HISTORY"));
                
                // Display detailed iterations if not too many
                if (Iterations.Count <= 20)
                {
                    foreach (var iteration in Iterations)
                    {
                        DisplayIterationDetails(iteration);
                    }
                }
                else
                {
                    Console.WriteLine(OutputFormatter.CreateBox(
                        "Detailed iteration output is available but was suppressed due to length.\n" +
                        "Consider increasing the verbosity level to see all iterations.",
                        "ITERATION DETAILS"
                    ));
                }
            }
            
            // Display final basis inverse if available and not too large
            if (BasisInverse != null && BasisInverse.GetLength(0) <= 10)
            {
                var basisInverseRows = new List<object[]>();
                for (int i = 0; i < BasisInverse.GetLength(0); i++)
                {
                    var row = new object[BasisInverse.GetLength(1) + 1];
                    row[0] = $"Row {i + 1}";
                    for (int j = 0; j < BasisInverse.GetLength(1); j++)
                    {
                        row[j + 1] = BasisInverse[i, j].ToString("F4");
                    }
                    basisInverseRows.Add(row);
                }
                
                var headers = new[] { "Basis Inverse" }
                    .Concat(Enumerable.Range(1, BasisInverse.GetLength(1)).Select(i => $"Col {i}"));
                    
                Console.WriteLine("FINAL BASIS INVERSE:");
                Console.WriteLine(OutputFormatter.CreateTable(headers.ToArray(), basisInverseRows));
            }
        }

        // Removed duplicate property definitions that were already in the ILinearProgramSolution implementation
        
        /// <summary>
        /// Displays detailed information about a single iteration
        /// </summary>
        private void DisplayIterationDetails(ISimplexIteration iteration)
        {
            var details = new StringBuilder();
            
            // Basic iteration info
            details.AppendLine(OutputFormatter.FormatKeyValue("Iteration", iteration.Iteration));
            details.AppendLine(OutputFormatter.FormatKeyValue("Status", iteration.Description));
            
            // Entering and leaving variables
            if (iteration.EnteringVariable >= 0 && VariableNames != null && 
                iteration.EnteringVariable < VariableNames.Length)
            {
                details.AppendLine(OutputFormatter.FormatKeyValue("Entering Variable", 
                    $"{VariableNames[iteration.EnteringVariable]} (x{iteration.EnteringVariable + 1})"));
            }
            
            if (iteration.LeavingVariable >= 0 && BasisVariables != null && 
                iteration.LeavingVariable < BasisVariables.Length)
            {
                details.AppendLine(OutputFormatter.FormatKeyValue("Leaving Variable", 
                    $"{BasisVariables[iteration.LeavingVariable]} (row {iteration.LeavingVariable + 1})"));
            }
            
            // Objective value
            details.AppendLine(OutputFormatter.FormatKeyValue("Objective Value", iteration.ObjectiveValue.ToString("F4")));
            
            // Pivot element if available
            if (iteration is TableauIteration tabIter && tabIter.PivotElement != 0)
            {
                details.AppendLine(OutputFormatter.FormatKeyValue("Pivot Element", tabIter.PivotElement.ToString("F4")));
            }
            
            Console.WriteLine(OutputFormatter.CreateBox(details.ToString(), $"ITERATION {iteration.Iteration}"));
            
            // Display basic variables and their values if available
            if (iteration is TableauIteration tabIter2 && tabIter2.BasicVariableValues != null)
            {
                var basicVars = tabIter2.BasisVariables
                    .Zip(tabIter2.BasicVariableValues, (name, value) => new object[] { name, value });
                
                Console.WriteLine("BASIC VARIABLES");
                Console.WriteLine(OutputFormatter.CreateTable(
                    new[] { "Basic Variable", "Value" },
                    basicVars
                ));
            }
        }

        // Removed duplicate property definitions that were already in the ILinearProgramSolution implementation

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
