using System;
using System.Linq;
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
                    var branchAndBoundSolver = new BranchAndBoundSolver();
                    
                    // Display canonical form first
                    Console.WriteLine("=== CANONICAL FORM ===");
                    branchAndBoundSolver.DisplayCanonicalForm(canonicalModel);
                    
                    Console.WriteLine("=== STARTING BRANCH AND BOUND ALGORITHM ===");
                    var bbSolution = branchAndBoundSolver.Solve(parsedModel);

                    // Display results
                    Console.WriteLine("\n=== BRANCH AND BOUND RESULTS ===");
                    Console.WriteLine($"Overall Status: {bbSolution.OverallStatus}");
                    Console.WriteLine($"Total Nodes Explored: {bbSolution.TotalNodesExplored}");
                    Console.WriteLine($"Total Nodes Fathomed: {bbSolution.TotalNodesFathomed}");

                    if (bbSolution.BestCandidate != null)
                    {
                        Console.WriteLine($"\n=== BEST INTEGER SOLUTION ===");
                        Console.WriteLine($"Objective Value: {bbSolution.BestObjectiveValue:F3}");
                        Console.WriteLine("Solution Vector:");
                        if (bbSolution.BestSolution != null)
                        {
                            for (int i = 0; i < bbSolution.BestSolution.Length; i++)
                            {
                                Console.WriteLine($"x{i + 1} = {bbSolution.BestSolution[i]:F3}");
                            }
                        }
                    }

                    // Display all nodes summary
                    Console.WriteLine("\n=== ALL NODES SUMMARY ===");
                    foreach (var node in bbSolution.AllNodes.OrderBy(n => n.Id))
                    {
                        string status = node.IsFathomed ? $"FATHOMED ({node.FathomReason})" : node.Solution?.Status ?? "Not Processed";
                        string objValue = node.Solution?.Status == "Optimal" ? $"{node.Solution.ObjectiveValue:F3}" : "N/A";
                        Console.WriteLine($"Node {node.Id}: {node.BranchingConstraint} → Status: {status}, Obj: {objValue}");
                    }
                }
                else
                {
                    // Regular simplex algorithms
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

                    if (hasIntegerVariables)
                    {
                        Console.WriteLine("\nNote: This is the LP relaxation solution. For integer solution, use Branch and Bound (option 3).");
                    }
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
