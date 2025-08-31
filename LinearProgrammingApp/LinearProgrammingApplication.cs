using System;
using System.Linq;
using LinearProgramming.Parsing;  // Needed for ParsedLinearProgrammingModel
using LinearProgramming.Algorithms;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.BranchAndBound;
using LinearProgramming.Algorithms.CuttingPlane;
using LinearProgramming.Algorithms.Knapsack;
using System.Collections.Generic;
// Use the nested CanonicalLinearProgrammingModel type from the correct namespace
using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;

namespace LinearProgrammingApp
{
    public class LinearProgrammingApplication
    {
        private readonly ParsedLinearProgrammingModel _model;
        private bool _isSolved = false;
        private object _solution = null;
        private string _algorithmUsed = "";

        public LinearProgrammingApplication(ParsedLinearProgrammingModel model)
        {
            _model = model ?? throw new ArgumentNullException(nameof(model));
        }

        public void Run()
        {
            bool exit = false;
            while (!exit)
            {
                Console.Clear();
                Console.WriteLine("=== Linear Programming Solver ===");
                Console.WriteLine("1. Solve Model");
                Console.WriteLine("2. View Model Information");
                Console.WriteLine("3. Sensitivity Analysis");
                Console.WriteLine("4. Exit");
                Console.Write("\nSelect an option: ");

                var key = Console.ReadKey();
                Console.WriteLine("\n");

                switch (key.KeyChar)
                {
                    case '1':
                        SolveModelMenu();
                        break;
                    case '2':
                        ViewModelInformation();
                        break;
                    case '3':
                        if (_isSolved)
                            ShowSensitivityAnalysis();
                        else
                            Console.WriteLine("Please solve the model first!");
                        break;
                    case '4':
                        exit = true;
                        break;
                    default:
                        Console.WriteLine("Invalid option. Please try again.");
                        break;
                }

                if (!exit)
                {
                    Console.WriteLine("\nPress any key to continue...");
                    Console.ReadKey();
                }
            }
        }

        private void SolveModelMenu()
        {
            bool hasIntegerVars = _model.Variables.Any(v => v.Type == VariableType.Integer || v.Type == VariableType.Binary);
            
            Console.WriteLine("\n=== Select Algorithm ===");
            Console.WriteLine("1. Primal Simplex");
            Console.WriteLine("2. Revised Primal Simplex");
            
            if (hasIntegerVars)
            {
                Console.WriteLine("3. Branch and Bound (for Integer Programming)");
                Console.WriteLine("4. Cutting Plane (for Integer Programming)");
            }
            
            Console.Write("\nSelect an algorithm: ");
            var key = Console.ReadKey();
            Console.WriteLine("\n");

            try
            {
                switch (key.KeyChar)
                {
                    case '1':
                        SolveWithPrimalSimplex();
                        break;
                    case '2':
                        SolveWithRevisedPrimalSimplex();
                        break;
                    case '3' when hasIntegerVars:
                        SolveWithBranchAndBound();
                        break;
                    case '4' when hasIntegerVars:
                        SolveWithCuttingPlane();
                        break;
                    default:
                        Console.WriteLine("Invalid selection or algorithm not applicable.");
                        return;
                }
                _isSolved = true;
                Console.WriteLine("\nSolution completed successfully!");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nError during solving: {ex.Message}");
            }
        }

        private void SolveWithPrimalSimplex()
        {
            var solver = new PrimalSimplexSolver();
            var canonicalModel = _model.ToCanonicalForm();
            _solution = solver.Solve(canonicalModel);
            _algorithmUsed = "Primal Simplex";
            DisplaySolution(_solution);
        }

        private void SolveWithRevisedPrimalSimplex()
        {
            var solver = new RevisedPrimalSimplexSolver();
            var canonicalModel = _model.ToCanonicalForm();
            _solution = solver.Solve(canonicalModel);
            _algorithmUsed = "Revised Primal Simplex";
            DisplaySolution(_solution);
        }

        private void SolveWithBranchAndBound()
        {
            var solver = new BranchAndBoundSolver();
            _solution = solver.Solve(_model);
            _algorithmUsed = "Branch and Bound";
            DisplaySolution(_solution);
        }

        private void SolveWithCuttingPlane()
        {
            var solver = new CuttingPlane();
            var canonicalModel = _model.ToCanonicalForm().ToTopLevelModel();
            _solution = solver.CuttingPlaneSolve(canonicalModel);
            _algorithmUsed = "Cutting Plane";
            DisplaySolution(_solution);
        }

        private void ViewModelInformation()
        {
            Console.WriteLine("\n=== Model Information ===");
            Console.WriteLine($"Objective: {_model.Objective.Optimization}");
            Console.WriteLine("\nVariables:");
            foreach (var variable in _model.Variables)
            {
                Console.WriteLine($"- {variable.Name}: {variable.Type}");
            }
            
            Console.WriteLine("\nConstraints:");
            for (int i = 0; i < _model.Constraints.Count; i++)
            {
                var constraint = _model.Constraints[i];
                Console.WriteLine($"{i + 1}. {constraint.Type} {constraint.RHS}");
            }
        }

        private void ShowSensitivityAnalysis()
        {
            if (!_isSolved || _solution == null)
            {
                Console.WriteLine("No solution available for sensitivity analysis.");
                return;
            }

            Console.WriteLine("\n=== Sensitivity Analysis ===");
            Console.WriteLine($"Algorithm used: {_algorithmUsed}");
            
            // Placeholder for actual sensitivity analysis
            // This would be implemented based on the specific solution type
            Console.WriteLine("\nSensitivity analysis results:");
            Console.WriteLine("- Objective value range: [TODO]");
            Console.WriteLine("- Shadow prices: [TODO]");
            Console.WriteLine("- Reduced costs: [TODO]");
        }

        private void DisplaySolution(object solution)
        {
            // This is a simplified display. In a real implementation, you would
            // handle different solution types appropriately
            Console.WriteLine("\n=== Solution ===");
            Console.WriteLine($"Algorithm: {_algorithmUsed}");
            
            if (solution is LinearProgramming.Algorithms.PrimalSimplex.LinearProgramSolution simplexSolution)
            {
                Console.WriteLine($"Status: {simplexSolution.Status}");
                Console.WriteLine($"Objective Value: {simplexSolution.ObjectiveValue:F4}");
                
                if (simplexSolution.SolutionVector != null)
                {
                    Console.WriteLine("\nSolution Vector:");
                    for (int i = 0; i < simplexSolution.SolutionVector.Length; i++)
                    {
                        Console.WriteLine($"x{i + 1} = {simplexSolution.SolutionVector[i]:F4}");
                    }
                }
            }
            // Add other solution type handlers as needed
        }
    }
}
