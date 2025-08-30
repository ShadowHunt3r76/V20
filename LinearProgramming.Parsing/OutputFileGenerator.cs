using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using LinearProgramming.Algorithms;

namespace LinearProgramming.Parsing
{
    /// <summary>
    /// Generates output files containing canonical form and tableau iterations for linear programming algorithms.
    /// All decimal values are rounded to three decimal places as required.
    /// </summary>
    public class OutputFileGenerator
    {
        /// <summary>
        /// Generates an output file for Primal Simplex algorithm results.
        /// </summary>
        /// <param name="parsedModel">The original parsed model</param>
        /// <param name="solution">The solution from Primal Simplex</param>
        /// <param name="outputPath">Output file path (optional, defaults to "primal_simplex_output.txt")</param>
        public static void GeneratePrimalSimplexOutput(ParsedLinearProgrammingModel parsedModel, LinearProgramSolution solution, string outputPath = "primal_simplex_output.txt")
        {
            using (StreamWriter writer = new StreamWriter(outputPath))
            {
                writer.WriteLine("=== PRIMAL SIMPLEX ALGORITHM OUTPUT ===");
                writer.WriteLine($"Generated on: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
                writer.WriteLine();

                // Write original problem
                WriteOriginalProblem(writer, parsedModel);

                // Write canonical form
                var canonicalModel = parsedModel.ToCanonicalForm();
                WriteCanonicalForm(writer, canonicalModel);

                // Write all tableau iterations
                WriteTableauIterations(writer, solution, "Primal Simplex");

                // Write final results
                WriteFinalResults(writer, solution);
            }
        }

        /// <summary>
        /// Generates an output file for Revised Primal Simplex algorithm results.
        /// </summary>
        /// <param name="parsedModel">The original parsed model</param>
        /// <param name="solution">The solution from Revised Primal Simplex</param>
        /// <param name="outputPath">Output file path (optional, defaults to "revised_simplex_output.txt")</param>
        public static void GenerateRevisedSimplexOutput(ParsedLinearProgrammingModel parsedModel, LinearProgramSolution solution, string outputPath = "revised_simplex_output.txt")
        {
            using (StreamWriter writer = new StreamWriter(outputPath))
            {
                writer.WriteLine("=== REVISED PRIMAL SIMPLEX ALGORITHM OUTPUT ===");
                writer.WriteLine($"Generated on: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
                writer.WriteLine();

                // Write original problem
                WriteOriginalProblem(writer, parsedModel);

                // Write canonical form
                var canonicalModel = parsedModel.ToCanonicalForm();
                WriteCanonicalForm(writer, canonicalModel);

                // Note: Revised simplex doesn't store full tableaux, so we note this
                writer.WriteLine("=== ALGORITHM NOTES ===");
                writer.WriteLine("Revised Primal Simplex uses matrix operations and doesn't maintain full tableaux.");
                writer.WriteLine("The algorithm works with basis matrices and reduced costs internally.");
                writer.WriteLine();

                // Write final results
                WriteFinalResults(writer, solution);
            }
        }

        /// <summary>
        /// Generates an output file for Branch and Bound algorithm results.
        /// </summary>
        /// <param name="parsedModel">The original parsed model</param>
        /// <param name="branchAndBoundSolution">The solution from Branch and Bound</param>
        /// <param name="outputPath">Output file path (optional, defaults to "branch_bound_output.txt")</param>
        public static void GenerateBranchAndBoundOutput(ParsedLinearProgrammingModel parsedModel, object branchAndBoundSolution, string outputPath = "branch_bound_output.txt")
        {
            using (StreamWriter writer = new StreamWriter(outputPath))
            {
                writer.WriteLine("=== BRANCH AND BOUND ALGORITHM OUTPUT ===");
                writer.WriteLine($"Generated on: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
                writer.WriteLine();

                // Write original problem
                WriteOriginalProblem(writer, parsedModel);

                // Write canonical form
                var canonicalModel = parsedModel.ToCanonicalForm();
                WriteCanonicalForm(writer, canonicalModel);

                // Write Branch and Bound specific results
                writer.WriteLine("=== BRANCH AND BOUND TREE ===");
                writer.WriteLine("Branch and Bound explores multiple nodes in a tree structure.");
                writer.WriteLine("Each node represents a sub-problem with additional constraints.");
                writer.WriteLine("Results show the exploration tree and final integer solution.");
                writer.WriteLine();

                // Note: The actual branch and bound solution structure would need to be accessed here
                // This is a placeholder that can be extended when the exact solution structure is available
                writer.WriteLine("=== FINAL INTEGER SOLUTION ===");
                writer.WriteLine("Branch and Bound results would be displayed here.");
            }
        }

        /// <summary>
        /// Generates a comprehensive output file that includes canonical form and tableau iterations.
        /// </summary>
        /// <param name="parsedModel">The original parsed model</param>
        /// <param name="solution">The algorithm solution</param>
        /// <param name="algorithmName">Name of the algorithm used</param>
        /// <param name="outputPath">Output file path</param>
        public static void GenerateComprehensiveOutput(ParsedLinearProgrammingModel parsedModel, LinearProgramSolution solution, string algorithmName, string outputPath)
        {
            using (StreamWriter writer = new StreamWriter(outputPath))
            {
                writer.WriteLine($"=== {algorithmName.ToUpper()} ALGORITHM OUTPUT ===");
                writer.WriteLine($"Generated on: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
                writer.WriteLine();

                // Write original problem
                WriteOriginalProblem(writer, parsedModel);

                // Write canonical form
                var canonicalModel = parsedModel.ToCanonicalForm();
                WriteCanonicalForm(writer, canonicalModel);

                // Write tableau iterations
                WriteTableauIterations(writer, solution, algorithmName);

                // Write final results
                WriteFinalResults(writer, solution);
            }
        }

        #region Private Helper Methods

        /// <summary>
        /// Writes the original problem formulation to the output file.
        /// </summary>
        private static void WriteOriginalProblem(StreamWriter writer, ParsedLinearProgrammingModel parsedModel)
        {
            writer.WriteLine("=== ORIGINAL PROBLEM ===");
            
            // Objective function
            string objType = parsedModel.Objective.Optimization == OptimizationType.Maximize ? "max" : "min";
            writer.Write($"{objType} ");
            
            for (int i = 0; i < parsedModel.Objective.Coefficients.Count; i++)
            {
                double coeff = parsedModel.Objective.Coefficients[i];
                string sign = coeff >= 0 ? "+" : "";
                writer.Write($"{sign}{coeff:F3} ");
            }
            writer.WriteLine();

            // Constraints
            for (int i = 0; i < parsedModel.Constraints.Count; i++)
            {
                var constraint = parsedModel.Constraints[i];
                
                for (int j = 0; j < constraint.Coefficients.Count; j++)
                {
                    double coeff = constraint.Coefficients[j];
                    string sign = coeff >= 0 ? "+" : "";
                    writer.Write($"{sign}{coeff:F3} ");
                }
                
                string opStr = constraint.Type == ConstraintType.LessThanOrEqual ? "<=" :
                              constraint.Type == ConstraintType.GreaterThanOrEqual ? ">=" : "=";
                writer.WriteLine($"{opStr} {constraint.RHS:F3}");
            }

            // Variable types
            for (int i = 0; i < parsedModel.Variables.Count; i++)
            {
                string varType = parsedModel.Variables[i].Type switch
                {
                    VariableType.NonNegative => "+",
                    VariableType.NonPositive => "-",
                    VariableType.Unrestricted => "urs",
                    VariableType.Integer => "int",
                    VariableType.Binary => "bin",
                    _ => "+"
                };
                writer.Write($"{varType} ");
            }
            writer.WriteLine();
            writer.WriteLine();
        }

        /// <summary>
        /// Writes the canonical form to the output file.
        /// </summary>
        private static void WriteCanonicalForm(StreamWriter writer, ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel canonicalModel)
        {
            writer.WriteLine("=== CANONICAL FORM ===");
            
            writer.WriteLine("Coefficient Matrix (A):");
            for (int i = 0; i < canonicalModel.CoefficientMatrix.Length; i++)
            {
                writer.Write("[ ");
                for (int j = 0; j < canonicalModel.CoefficientMatrix[i].Length; j++)
                {
                    writer.Write($"{canonicalModel.CoefficientMatrix[i][j]:F3,8} ");
                }
                writer.WriteLine("]");
            }
            writer.WriteLine();

            writer.WriteLine("RHS Vector (b):");
            writer.Write("[ ");
            for (int i = 0; i < canonicalModel.RHSVector.Length; i++)
            {
                writer.Write($"{canonicalModel.RHSVector[i]:F3,8} ");
            }
            writer.WriteLine("]");
            writer.WriteLine();

            writer.WriteLine("Objective Coefficients (c):");
            writer.Write("[ ");
            for (int i = 0; i < canonicalModel.ObjectiveCoefficients.Length; i++)
            {
                writer.Write($"{canonicalModel.ObjectiveCoefficients[i]:F3,8} ");
            }
            writer.WriteLine("]");
            writer.WriteLine();

            writer.WriteLine("Variable Types:");
            for (int i = 0; i < canonicalModel.VariableTypes.Length; i++)
            {
                string varType = canonicalModel.VariableTypes[i] switch
                {
                    VariableType.NonNegative => "NonNegative",
                    VariableType.NonPositive => "NonPositive", 
                    VariableType.Unrestricted => "Unrestricted",
                    VariableType.Integer => "Integer",
                    VariableType.Binary => "Binary",
                    _ => "NonNegative"
                };
                writer.WriteLine($"x{i + 1}: {varType}");
            }
            writer.WriteLine();
        }

        /// <summary>
        /// Writes all tableau iterations to the output file.
        /// </summary>
        private static void WriteTableauIterations(StreamWriter writer, LinearProgramSolution solution, string algorithmName)
        {
            writer.WriteLine($"=== {algorithmName.ToUpper()} TABLEAU ITERATIONS ===");
            
            if (solution.InitialTable != null)
            {
                writer.WriteLine("Initial Tableau:");
                WriteTableau(writer, solution.InitialTable);
                writer.WriteLine();
            }

            if (solution.TableHistory != null && solution.TableHistory.Count > 0)
            {
                for (int i = 0; i < solution.TableHistory.Count; i++)
                {
                    writer.WriteLine($"Iteration {i + 1}:");
                    WriteTableau(writer, solution.TableHistory[i]);
                    writer.WriteLine();
                }
            }

            if (solution.OptimalTable != null)
            {
                writer.WriteLine("Final Tableau:");
                WriteTableau(writer, solution.OptimalTable);
                writer.WriteLine();
            }
        }

        /// <summary>
        /// Writes a single tableau to the output file with proper formatting.
        /// </summary>
        private static void WriteTableau(StreamWriter writer, double[,] tableau)
        {
            if (tableau == null) return;

            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);

            // Write tableau header
            writer.Write("     ");
            for (int j = 0; j < cols - 1; j++)
            {
                writer.Write($"{"x" + (j + 1),8} ");
            }
            writer.WriteLine($"{"RHS",8}");

            // Write separator line
            writer.Write("     ");
            for (int j = 0; j < cols; j++)
            {
                writer.Write("---------");
            }
            writer.WriteLine();

            // Write tableau rows
            for (int i = 0; i < rows; i++)
            {
                if (i == 0)
                    writer.Write("Z  | ");
                else
                    writer.Write($"R{i} | ");

                for (int j = 0; j < cols; j++)
                {
                    writer.Write($"{tableau[i, j]:F3,8} ");
                }
                writer.WriteLine();
            }
        }

        /// <summary>
        /// Writes the final solution results to the output file.
        /// </summary>
        private static void WriteFinalResults(StreamWriter writer, LinearProgramSolution solution)
        {
            writer.WriteLine("=== FINAL RESULTS ===");
            writer.WriteLine($"Status: {solution.Status}");
            writer.WriteLine($"Objective Value: {solution.ObjectiveValue:F3}");
            writer.WriteLine();

            if (solution.SolutionVector != null)
            {
                writer.WriteLine("Solution Vector:");
                for (int i = 0; i < solution.SolutionVector.Length; i++)
                {
                    writer.WriteLine($"x{i + 1} = {solution.SolutionVector[i]:F3}");
                }
            }
            writer.WriteLine();
            writer.WriteLine("=== END OF OUTPUT ===");
        }

        #endregion
    }
}