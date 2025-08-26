using System;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms;

namespace LinearProgrammingApp
{
    class Program
    {
        static void Main(string[] args)
        {
            // Example: Read input file path from args
            if (args.Length == 0)
            {
                Console.WriteLine("Usage: PrimalSolver <inputfile>");
                return;
            }
            string inputFile = args[0];
            try
            {
                var parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
                var parsedModel = parser.ParseFromFile(inputFile);
                var canonicalModel = parsedModel.ToCanonicalForm();

                // Choose algorithm (Primal or Revised)
                Console.WriteLine("Choose algorithm: 1) Primal Simplex  2) Revised Primal Simplex");
                var key = Console.ReadKey();
                Console.WriteLine();
                LinearProgramSolution solution = null;
                if (key.KeyChar == '2')
                {
                    var revisedSolver = new RevisedPrimalSimplexSolver();
                    solution = revisedSolver.Solve(canonicalModel);
                }
                else
                {
                    var primalSolver = new PrimalSimplexSolver();
                    solution = primalSolver.Solve(canonicalModel);
                }

                // Output results
                Console.WriteLine($"Status: {solution.Status}");
                if (solution.SolutionVector != null)
                {
                    Console.WriteLine("Solution vector:");
                    for (int i = 0; i < solution.SolutionVector.Length; i++)
                        Console.WriteLine($"x{i + 1} = {solution.SolutionVector[i]:F4}");
                }
                Console.WriteLine($"Objective value: {solution.ObjectiveValue:F4}");
            }
            catch (ParsedLinearProgrammingModel.LinearProgrammingParseException ex)
            {
                Console.WriteLine($"Input error: {ex.Message}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Unexpected error: {ex.Message}");
            }
        }
    }
}
