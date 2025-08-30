using System;
using System.Linq;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms;
using LinearProgramming.Algorithms.BranchAndBound;
using LinearProgramming.Algorithms.PrimalSimplex;

namespace LinearProgrammingApp
{
    class Program
    {
        static void Main(string[] args)
        {
            // Example: Read input file path from args
            if (args.Length == 0)
            {
                Console.WriteLine("Usage: LinearProgrammingApp <inputfile>");
                return;
            }
            
            try
            {
                string inputFile = args[0];
                var parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
                var parsedModel = parser.ParseFromFile(inputFile);
                var canonicalModel = parsedModel.ToCanonicalForm();

                //Call menu
                var solverMenu = new SolverMenu();
                solverMenu.RunMenu(parsedModel);

                // Check if the problem has integer/binary variables
                bool hasIntegerVariables = parsedModel.Variables.Any(v => v.Type == VariableType.Integer || v.Type == VariableType.Binary);

                if (hasIntegerVariables)
                {
                    Console.WriteLine("Integer Programming Problem Detected!");
                    Console.WriteLine("Choose algorithm:");
                    Console.WriteLine("1) Primal Simplex (LP Relaxation)");
                    Console.WriteLine("2) Revised Primal Simplex (LP Relaxation)"); 
                    Console.WriteLine("3) Branch and Bound Simplex Algorithm");
                }
                else
                {
                    Console.WriteLine("Linear Programming Problem Detected.");
                    Console.WriteLine("Choose algorithm:");
                    Console.WriteLine("1) Primal Simplex");
                    Console.WriteLine("2) Revised Primal Simplex");
                }

                var key = Console.ReadKey();
                Console.WriteLine();
                Console.WriteLine();

                if (hasIntegerVariables && key.KeyChar == '3')
                {
                    // Branch and Bound Algorithm
                    var solver = new SolverMenu();
                    solver.RunBranchAndBound(parsedModel);
                    return;
                }
                
                // Regular simplex algorithms
                LinearProgramSolution solution = null;
                
                if (key.KeyChar == '2')
                {
                    var revisedSolver = new RevisedPrimalSimplexSolver();
                    solution = revisedSolver.Solve(canonicalModel);
                    
                    // Generate output file
                    OutputFileGenerator.GenerateRevisedSimplexOutput(parsedModel, solution, "revised_simplex_main_output.txt");
                    Console.WriteLine("\nOutput file generated: revised_simplex_main_output.txt");
                }
                else
                {
                    var primalSolver = new PrimalSimplexSolver();
                    solution = primalSolver.Solve(canonicalModel);
                    
                    // Generate output file
                    OutputFileGenerator.GeneratePrimalSimplexOutput(parsedModel, solution, "primal_simplex_main_output.txt");
                    Console.WriteLine("\nOutput file generated: primal_simplex_main_output.txt");
                }

                // Output results to console
                Console.WriteLine($"Status: {solution.Status}");
                if (solution.SolutionVector != null)
                {
                    Console.WriteLine("Solution vector:");
                    for (int i = 0; i < solution.SolutionVector.Length; i++)
                        Console.WriteLine($"x{i + 1} = {solution.SolutionVector[i]:F3}");
                }
                Console.WriteLine($"Objective value: {solution.ObjectiveValue:F3}");

                if (hasIntegerVariables)
                {
                    Console.WriteLine("\nNote: This is the LP relaxation solution. For integer solution, use Branch and Bound (option 3).");
                }
            }
            catch (ParsedLinearProgrammingModel.LinearProgrammingParseException ex)
            {
                Console.WriteLine($"Input error: {ex.Message}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Unexpected error: {ex.Message}");
                Console.WriteLine($"Stack trace: {ex.StackTrace}");
            }
        }
    }
}
