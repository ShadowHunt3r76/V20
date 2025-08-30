using System;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms;

namespace LinearProgrammingApp
{
    public class SolverMenu

        //Interactive menu for application to select chosen solver
    {
        public void RunMenu(ParsedLinearProgrammingModel parsedModel)
        {
            bool exit = false;

            while (!exit)
            {
                Console.WriteLine("\n=== Solver Menu ===");
                Console.WriteLine("1) Primal Simplex");
                Console.WriteLine("2) Revised Primal Simplex");
                Console.WriteLine("3) Cutting Plane");
                Console.WriteLine("4) Branch and Bound");
                Console.WriteLine("5) Knapsack");
                Console.WriteLine("0) Exit");
                Console.Write("Choose an option: ");

                string choice = Console.ReadLine();
                Console.WriteLine();

                switch (choice)
                {
                    case "1":
                        RunPrimalSimplex(parsedModel);
                        break;

                    case "2":
                        RunRevisedPrimalSimplex(parsedModel);
                        break;

                    case "3":
                        RunCuttingPlane(parsedModel);
                        break;

                    case "4":
                        RunBranchAndBound(parsedModel);
                        break;

                    case "5":
                        RunKnapsack();
                        break;

                    case "0":
                        exit = true;
                        Console.WriteLine("Exiting solver menu...");
                        break;

                    default:
                        Console.WriteLine("Invalid choice, please try again.");
                        break;
                }
            }
        }

        private void RunPrimalSimplex(ParsedLinearProgrammingModel parsedModel)
        {
            var canonicalModel = parsedModel.ToCanonicalForm();
            var solver = new PrimalSimplexSolver();
            var solution = solver.Solve(canonicalModel);

            Console.WriteLine("\n--- Primal Simplex Result ---");
            Console.WriteLine($"Status: {solution.Status}");
            if (solution.SolutionVector != null)
            {
                for (int i = 0; i < solution.SolutionVector.Length; i++)
                {
                    Console.WriteLine($"x{i + 1} = {solution.SolutionVector[i]:F3}");
                }
            }
            Console.WriteLine($"Objective Value = {solution.ObjectiveValue:F3}");

            // Generate output file
            OutputFileGenerator.GeneratePrimalSimplexOutput(parsedModel, solution);
            Console.WriteLine("\nOutput file generated: primal_simplex_output.txt");
        }

        private void RunRevisedPrimalSimplex(ParsedLinearProgrammingModel parsedModel)
        {
            var canonicalModel = parsedModel.ToCanonicalForm();
            var solver = new RevisedPrimalSimplexSolver();
            var solution = solver.Solve(canonicalModel);

            Console.WriteLine("\n--- Revised Primal Simplex Result ---");
            Console.WriteLine($"Status: {solution.Status}");
            if (solution.SolutionVector != null)
            {
                for (int i = 0; i < solution.SolutionVector.Length; i++)
                {
                    Console.WriteLine($"x{i + 1} = {solution.SolutionVector[i]:F3}");
                }
            }
            Console.WriteLine($"Objective Value = {solution.ObjectiveValue:F3}");

            // Generate output file
            OutputFileGenerator.GenerateRevisedSimplexOutput(parsedModel, solution);
            Console.WriteLine("\nOutput file generated: revised_simplex_output.txt");
        }

        private void RunCuttingPlane(ParsedLinearProgrammingModel parsedModel)
        {
            var canonicalModel = parsedModel.ToCanonicalForm();
            var solver = new CuttingPlaneAlgorithm.CuttingPlane();
            var solution = solver.CuttingPlaneSolve(canonicalModel);

            Console.WriteLine("\n=== CUTTING PLANE RESULTS ===");
            Console.WriteLine($"Status: {solution.Status}");
            Console.WriteLine($"Objective value: {solution.ObjectiveValue:F3}");

            if (solution.SolutionVector != null)
            {
                Console.WriteLine("Solution vector:");
                for (int i = 0; i < solution.SolutionVector.Length; i++)
                {
                    Console.WriteLine($"x{i + 1} = {solution.SolutionVector[i]:F3}");
                }
            }

            Console.WriteLine("\nTableau iterations were written to cuttingplane_output.txt");
            
            // Also generate using our comprehensive output format
            OutputFileGenerator.GenerateComprehensiveOutput(parsedModel, solution, "Cutting Plane", "cutting_plane_comprehensive_output.txt");
            Console.WriteLine("Comprehensive output file generated: cutting_plane_comprehensive_output.txt");
        }

        private void RunBranchAndBound(ParsedLinearProgrammingModel parsedModel)
        {
            var solver = new BranchAndBoundSolver();
            var solution = solver.Solve(parsedModel);

            Console.WriteLine("\n--- Branch and Bound Result ---");
            Console.WriteLine($"Overall Status: {solution.OverallStatus}");
            if (solution.BestSolution != null)
            {
                for (int i = 0; i < solution.BestSolution.Length; i++)
                {
                    Console.WriteLine($"x{i + 1} = {solution.BestSolution[i]:F3}");
                }
            }
            Console.WriteLine($"Best Objective Value = {solution.BestObjectiveValue:F3}");
            Console.WriteLine($"Total Nodes Explored: {solution.TotalNodesExplored}");
            Console.WriteLine($"Total Nodes Fathomed: {solution.TotalNodesFathomed}");

            // Generate output file for Branch and Bound
            OutputFileGenerator.GenerateBranchAndBoundOutput(parsedModel, solution);
            Console.WriteLine("\nOutput file generated: branch_bound_output.txt");
        }

        private void RunKnapsack()
        {
            Console.WriteLine("\n--- Knapsack Problem ---");
            Console.Write("Enter number of items: ");
            int n = int.Parse(Console.ReadLine());

            double[] weights = new double[n];
            double[] values = new double[n];

            for (int i = 0; i < n; i++)
            {
                Console.Write($"Enter weight of item {i + 1}: ");
                weights[i] = double.Parse(Console.ReadLine());

                Console.Write($"Enter value of item {i + 1}: ");
                values[i] = double.Parse(Console.ReadLine());
            }

            Console.Write("Enter capacity of knapsack: ");
            double capacity = double.Parse(Console.ReadLine());

            var knapsackSolver = new LPRProject.Knapsack();
            double maxValue = knapsackSolver.SolveKnapsack(weights, values, capacity);

            Console.WriteLine($"Maximum value in knapsack = {maxValue:F3}");
        }
    }
}
