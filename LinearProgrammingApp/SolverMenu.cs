using System;
using System.Linq;
using System.IO;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.CuttingPlane;
using LinearProgramming.Algorithms.BranchAndBound;
using LinearProgramming.Algorithms.Knapsack;

namespace LinearProgrammingApp
{
    /// <summary>
    /// Handles the solver selection and execution for different LP and IP problems
    /// </summary>
    public static class SolverMenu
    {
        /// <summary>
        /// Main menu for the Linear Programming Solver
        /// </summary>
        /// <param name="model">The parsed linear programming model</param>
        /// <param name="outputPath">Path to save the output file</param>
        public static void ShowMainMenu(ParsedLinearProgrammingModel model, string outputPath)
        {
            while (true)
            {
                Console.Clear();
                Console.WriteLine("=== Linear Programming Solver ===");
                Console.WriteLine("1. View Model Information");
                Console.WriteLine("2. Solve Model");
                Console.WriteLine("3. Sensitivity Analysis");
                Console.WriteLine("4. Exit");
                Console.Write("\nSelect an option: ");

                switch (Console.ReadLine())
                {
                    case "1":
                        DisplayModelInfo(model);
                        break;
                        
                    case "2":
                        ShowSolverMenu(model, outputPath);
                        break;
                        
                    case "3":
                        ShowSensitivityAnalysis(model);
                        break;
                        
                    case "4":
                        Console.WriteLine("\nExiting...");
                        return;
                        
                    default:
                        Console.WriteLine("\nInvalid option. Please try again.");
                        break;
                }

                if (Console.ReadKey(true).Key == ConsoleKey.Escape)
                    break;
            }
        }

        private static void ShowSolverMenu(ParsedLinearProgrammingModel model, string outputPath)
        {
            bool hasIntegerVars = model.Variables.Any(v => v.Type == VariableType.Integer || v.Type == VariableType.Binary);
            
            while (true)
            {
                Console.Clear();
                Console.WriteLine("=== Select Solver ===");
                Console.WriteLine("1. Primal Simplex (LP)");
                Console.WriteLine("2. Revised Primal Simplex (LP)");
                
                if (hasIntegerVars)
                {
                    Console.WriteLine("3. Branch and Bound (IP)");
                    Console.WriteLine("4. Cutting Plane (IP)");
                    Console.WriteLine("5. Knapsack Solver (Binary IP)");
                }
                else
                {
                    Console.WriteLine("3. Knapsack Solver");
                }
                
                Console.WriteLine("0. Back to Main Menu");
                Console.Write("\nSelect an option: ");

                string choice = Console.ReadLine();
                Console.WriteLine();

                try
                {
                    switch (choice)
                    {
                        case "1":
                            RunPrimalSimplex(model, outputPath);
                            break;
                            
                        case "2":
                            RunRevisedPrimalSimplex(model, outputPath);
                            break;
                            
                        case "3" when hasIntegerVars:
                            RunBranchAndBound(model, outputPath);
                            break;
                            
                        case "4" when hasIntegerVars:
                            RunCuttingPlane(model, outputPath);
                            break;
                            
                        case "5" when hasIntegerVars:
                        case "3" when !hasIntegerVars:
                            RunKnapsackSolver(model, outputPath);
                            break;
                            
                        case "0":
                            return;
                            
                        default:
                            Console.WriteLine("Invalid option. Please try again.");
                            break;
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"\nError: {ex.Message}");
                }

                Console.WriteLine("\nPress any key to continue...");
                Console.ReadKey();
            }
        }

        private static void RunPrimalSimplex(ParsedLinearProgrammingModel model, string outputPath)
        {
            Console.WriteLine("\n=== Primal Simplex Algorithm ===");
            var canonicalModel = model.ToCanonicalForm();
            var solver = new PrimalSimplexSolver();
            
            Console.WriteLine("Solving...");
            var solution = solver.Solve(canonicalModel);

            // Display results
            Console.WriteLine("\n=== Solution ===");
            Console.WriteLine($"Status: {solution.Status}");
            Console.WriteLine($"Objective Value: {solution.ObjectiveValue:F6}");
            
            if (solution.SolutionVector != null)
            {
                Console.WriteLine("\nVariable Values:");
                for (int i = 0; i < solution.SolutionVector.Length; i++)
                {
                    string varName = i < model.Variables.Count ? model.Variables[i].Name : $"s{i - model.Variables.Count + 1}";
                    Console.WriteLine($"{varName} = {solution.SolutionVector[i]:F6}");
                }
            }

            // Save results
            string solutionPath = Path.Combine(Path.GetDirectoryName(outputPath), $"primal_simplex_{Path.GetFileName(outputPath)}");
            OutputFileGenerator.GenerateSimplexOutput(model, solution, "Primal Simplex", solutionPath);
            Console.WriteLine($"\nSolution saved to: {Path.GetFullPath(solutionPath)}");
        }

        private static void RunRevisedPrimalSimplex(ParsedLinearProgrammingModel model, string outputPath)
        {
            Console.WriteLine("\n=== Revised Primal Simplex Algorithm ===");
            var canonicalModel = model.ToCanonicalForm();
            var solver = new RevisedPrimalSimplexSolver();
            
            Console.WriteLine("Solving...");
            var solution = solver.Solve(canonicalModel);

            // Display results
            Console.WriteLine("\n=== Solution ===");
            Console.WriteLine($"Status: {solution.Status}");
            Console.WriteLine($"Objective Value: {solution.ObjectiveValue:F6}");
            
            if (solution.SolutionVector != null)
            {
                Console.WriteLine("\nVariable Values:");
                for (int i = 0; i < solution.SolutionVector.Length; i++)
                {
                    string varName = i < model.Variables.Count ? model.Variables[i].Name : $"s{i - model.Variables.Count + 1}";
                    Console.WriteLine($"{varName} = {solution.SolutionVector[i]:F6}");
                }
            }

            // Save results
            string solutionPath = Path.Combine(Path.GetDirectoryName(outputPath), $"revised_primal_simplex_{Path.GetFileName(outputPath)}");
            OutputFileGenerator.GenerateSimplexOutput(model, solution, "Revised Primal Simplex", solutionPath);
            Console.WriteLine($"\nSolution saved to: {Path.GetFullPath(solutionPath)}");
        }

        private static void RunCuttingPlane(ParsedLinearProgrammingModel model, string outputPath)
        {
            Console.WriteLine("\n=== Cutting Plane Algorithm (Integer Programming) ===");
            var canonicalModel = model.ToCanonicalForm();
            var solver = new CuttingPlane();
            
            Console.WriteLine("Solving with Cutting Plane method...");
            var solution = solver.CuttingPlaneSolve(canonicalModel.ToTopLevelModel());

            // Display results
            Console.WriteLine("\n=== Integer Solution ===");
            Console.WriteLine($"Status: {solution.Status}");
            Console.WriteLine($"Objective Value: {solution.ObjectiveValue:F0}");
            
            if (solution.SolutionVector != null)
            {
                Console.WriteLine("\nVariable Values:");
                for (int i = 0; i < solution.SolutionVector.Length; i++)
                {
                    string varName = i < model.Variables.Count ? model.Variables[i].Name : $"s{i - model.Variables.Count + 1}";
                    Console.WriteLine($"{varName} = {solution.SolutionVector[i]:F0}");
                }
            }

            // Save results
            string solutionPath = Path.Combine(Path.GetDirectoryName(outputPath), $"cutting_plane_{Path.GetFileName(outputPath)}");
            OutputFileGenerator.GenerateIntegerSolutionOutput(model, solution, "Cutting Plane", solutionPath);
            Console.WriteLine($"\nSolution saved to: {Path.GetFullPath(solutionPath)}");
        }

        private static void RunBranchAndBound(ParsedLinearProgrammingModel model, string outputPath)
        {
            Console.WriteLine("\n=== Branch and Bound Algorithm (Integer Programming) ===");
            var solver = new BranchAndBoundSolver();
            
            Console.WriteLine("Solving with Branch and Bound...");
            var solution = solver.Solve(model);

            // Display results
            Console.WriteLine("\n=== Integer Solution ===");
            Console.WriteLine($"Status: {solution.OverallStatus}");
            Console.WriteLine($"Objective Value: {solution.BestObjectiveValue:F0}");
            Console.WriteLine($"Nodes Explored: {solution.TotalNodesExplored}");
            Console.WriteLine($"Nodes Fathomed: {solution.TotalNodesFathomed}");
            
            if (solution.BestSolution != null)
            {
                Console.WriteLine("\nVariable Values:");
                for (int i = 0; i < solution.BestSolution.Length; i++)
                {
                    string varName = i < model.Variables.Count ? model.Variables[i].Name : $"s{i - model.Variables.Count + 1}";
                    Console.WriteLine($"{varName} = {solution.BestSolution[i]:F0}");
                }
            }

            // Save results
            string solutionPath = Path.Combine(Path.GetDirectoryName(outputPath), $"branch_and_bound_{Path.GetFileName(outputPath)}");
            OutputFileGenerator.GenerateBranchAndBoundOutput(model, solution, solutionPath);
            Console.WriteLine($"\nSolution saved to: {Path.GetFullPath(solutionPath)}");
        }

        private static void RunKnapsackSolver(ParsedLinearProgrammingModel model, string outputPath)
        {
            Console.WriteLine("\n=== Knapsack Problem Solver ===");
            
            try 
            {
                // Extract values from the objective function (maximization)
                double[] values = model.Objective.Coefficients.ToArray();
                int n = values.Length;
                
                // Extract weights from the first constraint
                double[] weights = model.Constraints[0].Coefficients.ToArray();
                double capacity = model.Constraints[0].RHS;
                
                // Display the problem being solved
                Console.WriteLine("\nSolving knapsack problem with:");
                Console.WriteLine($"- Items: {n}");
                Console.WriteLine($"- Values: {string.Join(", ", values)}");
                Console.WriteLine($"- Weights: {string.Join(", ", weights)}");
                Console.WriteLine($"- Capacity: {capacity}");

                // Create knapsack problem with detailed logging
                var knapsack = new LinearProgramming.Algorithms.Knapsack.KnapsackSolver(weights, values, capacity);

                // Enable detailed logging for branch and bound
                Console.WriteLine("\n=== Branch and Bound Process ===");
                var (maxValue, includedItems) = knapsack.Solve(enableLogging: true);

                // Display results
                Console.WriteLine("\n=== Optimal Solution ===");
                Console.WriteLine($"Maximum value: {maxValue:F2}");
                
                double totalWeight = 0;
                Console.WriteLine("\nSelected Items:");
                for (int i = 0; i < includedItems.Length; i++)
                {
                    if (includedItems[i])
                    {
                        totalWeight += weights[i];
                        Console.WriteLine($"  Item {i + 1}: Value = {values[i]:F2}, Weight = {weights[i]:F2}");
                    }
                }
                Console.WriteLine($"\nTotal Weight: {totalWeight:F2}/{capacity}");
                
                // Save results
                string solutionPath = Path.Combine(Path.GetDirectoryName(outputPath), 
                    $"knapsack_{Path.GetFileName(outputPath)}");
                OutputFileGenerator.GenerateKnapsackOutput(weights, values, capacity, (maxValue, includedItems), solutionPath);
                Console.WriteLine($"\nSolution saved to: {Path.GetFullPath(solutionPath)}");
                
                // Ask if user wants to perform sensitivity analysis
                Console.Write("\nPerform sensitivity analysis? (y/n): ");
                if (Console.ReadLine().Trim().ToLower() == "y")
                {
                    PerformKnapsackSensitivityAnalysis(weights, values, capacity, maxValue, includedItems, solutionPath);
                }
                
                // Display the complete search tree
                Console.WriteLine("\n=== Complete Branch and Bound Search Tree ===");
                DisplayKnapsackSearchTree(knapsack);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nError solving knapsack problem: {ex.Message}");
                Console.WriteLine(ex.StackTrace);
            }
            
            Console.WriteLine("\nPress any key to continue...");
            Console.ReadKey();
        }
        
        private static void DisplayKnapsackSearchTree(LinearProgramming.Algorithms.Knapsack.KnapsackSolver knapsack)
        {
            try
            {
                var searchTree = knapsack.SearchTree;
                if (searchTree == null || searchTree.Count == 0)
                {
                    Console.WriteLine("No search tree data available.");
                    return;
                }
                
                Console.WriteLine("\nSearch Tree Visualization:");
                Console.WriteLine("=========================\n");
                
                // Group nodes by level for better visualization
                var levels = searchTree.GroupBy(n => n.Level).OrderBy(g => g.Key);
                
                foreach (var level in levels)
                {
                    Console.WriteLine($"Level {level.Key + 1}:");
                    foreach (var node in level)
                    {
                        string nodeInfo = $"  {node.Path} [Value: {node.Value}, Bound: {node.Bound:F2}, ";
                        nodeInfo += $"Weight: {node.Weight}] ";
                        
                        // Highlight the best solution
                        if (node == knapsack.BestNode)
                        {
                            Console.ForegroundColor = ConsoleColor.Green;
                            nodeInfo += "(Best Solution)";
                        }
                        
                        Console.WriteLine(nodeInfo);
                        Console.ResetColor();
                        
                        // Show included items
                        var included = node.Included.Select((inc, i) => inc ? i + 1 : -1)
                                                 .Where(i => i > 0)
                                                 .ToList();
                        
                        if (included.Any())
                        {
                            Console.WriteLine($"    Included items: {string.Join(", ", included)}");
                        }
                    }
                    Console.WriteLine();
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error displaying search tree: {ex.Message}");
            }
        }

        private static void DisplayModelInfo(ParsedLinearProgrammingModel model)
        {
            Console.Clear();
            Console.WriteLine("=== Model Information ===");
            
            // Basic info
            Console.WriteLine($"\nObjective: {(model.Objective.Optimization == OptimizationType.Maximize ? "Maximize" : "Minimize")} {string.Join(" + ", model.Objective.Coefficients.Select((c, i) => $"{c}x{i+1}"))}");
            
            // Variables
            Console.WriteLine("\nVariables:");
            for (int i = 0; i < model.Variables.Count; i++)
            {
                var variable = model.Variables[i];
                string type = variable.Type switch
                {
                    VariableType.NonNegative => "Non-negative",
                    VariableType.NonPositive => "Non-positive",
                    VariableType.Binary => "Binary (0 or 1)",
                    VariableType.Integer => "Integer",
                    VariableType.Unrestricted => "Unrestricted",
                    _ => "Unknown"
                };
                Console.WriteLine($"  x{i+1} : {type}");
            }
            
            // Constraints
            Console.WriteLine("\nConstraints:");
            foreach (var constraint in model.Constraints)
            {
                string lhs = string.Join(" + ", constraint.Coefficients.Select((c, i) => $"{c}x{i+1}"));
                string relation = constraint.Type == ConstraintType.LessThanOrEqual ? "<=" : ">=";
                Console.WriteLine($"  {lhs} {relation} {constraint.RHS}");
            }
            
            // Variable Bounds
            Console.WriteLine("\nVariable Bounds:");
            for (int i = 0; i < model.Variables.Count; i++)
            {
                var variable = model.Variables[i];
                string bounds = variable.Type switch
                {
                    VariableType.NonNegative => "0 ≤ x{i+1} < ∞",
                    VariableType.NonPositive => "-∞ < x{i+1} ≤ 0",
                    VariableType.Binary => "x{i+1} ∈ {0, 1}",
                    VariableType.Integer => "x{i+1} ∈ ℤ",
                    VariableType.Unrestricted => "-∞ < x{i+1} < ∞",
                    _ => "Unknown bounds"
                };
                Console.WriteLine($"  {bounds}");
            }
            
            Console.WriteLine("\nPress any key to continue...");
        }
        
        private static string GetRelationSymbol(ConstraintType type)
        {
            return type switch
            {
                ConstraintType.LessThanOrEqual => "<=",
                ConstraintType.GreaterThanOrEqual => ">=",
                ConstraintType.Equal => "=",
                _ => "?"
            };
        }
        
        private static string GetBoundType(double? lower, double? upper)
        {
            if (lower.HasValue && upper.HasValue)
                return $"in [{lower}, {upper}]";
            if (lower.HasValue)
                return $">= {lower}";
            if (upper.HasValue)
                return $"<= {upper}";
            return "free";
        }
        
        private static void ShowSensitivityAnalysis(ParsedLinearProgrammingModel model)
        {
            // For now, we'll assume no integer variables since we don't have that information
            bool hasIntegerVars = false;
            
            while (true)
            {
                Console.Clear();
                Console.WriteLine("=== Sensitivity Analysis ===");
                Console.WriteLine("Select the algorithm for sensitivity analysis:");
                Console.WriteLine("1. Primal Simplex");
                Console.WriteLine("2. Revised Primal Simplex");
                
                if (hasIntegerVars)
                {
                    Console.WriteLine("3. Branch and Bound");
                    Console.WriteLine("4. Cutting Plane");
                    Console.WriteLine("5. Knapsack Problem");
                }
                else
                {
                    Console.WriteLine("3. Knapsack Problem");
                }
                
                Console.WriteLine("0. Back to Main Menu");
                Console.Write("\nSelect an option: ");

                string choice = Console.ReadLine();
                Console.WriteLine();

                try
                {
                    switch (choice)
                    {
                        case "1":
                            RunPrimalSimplexSensitivity(model);
                            break;
                            
                        case "2":
                            RunRevisedPrimalSimplexSensitivity(model);
                            break;
                            
                        case "3" when hasIntegerVars:
                            RunBranchAndBoundSensitivity(model);
                            break;
                            
                        case "4" when hasIntegerVars:
                            RunCuttingPlaneSensitivity(model);
                            break;
                            
                        case "5" when hasIntegerVars:
                        case "3" when !hasIntegerVars:
                            RunKnapsackSensitivity();
                            break;
                            
                        case "0":
                            return;
                            
                        default:
                            Console.WriteLine("Invalid option. Please try again.");
                            break;
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"\nError during sensitivity analysis: {ex.Message}");
                    Console.WriteLine(ex.StackTrace);
                }

                Console.WriteLine("\nPress any key to continue...");
                Console.ReadKey();
            }
        }

        private static void RunPrimalSimplexSensitivity(ParsedLinearProgrammingModel model)
        {
            Console.WriteLine("\n=== Primal Simplex Sensitivity Analysis ===");
            var canonicalModel = model.ToCanonicalForm();
            var solver = new PrimalSimplexSolver();
            
            Console.WriteLine("Solving with Primal Simplex...");
            var solution = solver.Solve(canonicalModel);
            
            if (!string.Equals(solution.Status, "Optimal", StringComparison.OrdinalIgnoreCase))
            {
                Console.WriteLine("Cannot perform sensitivity analysis - the problem does not have an optimal solution.");
                return;
            }
            
            // Get constraint types and RHS from the canonical model
            var constraintTypes = canonicalModel.ConstraintTypes;
            var rhsCoefficients = canonicalModel.RHSVector;
            
            // Initialize sensitivity analysis with required parameters
            var sensitivity = new SimplexSensitivityAnalysis(
                finalTableau: solution.OptimalTable,
                basisIndices: solution.BasisIndices,
                variableNames: solution.VariableNames,
                objectiveCoefficients: canonicalModel.ObjectiveCoefficients,
                rhsCoefficients: rhsCoefficients,
                constraintTypes: constraintTypes
            );
            
            // Perform sensitivity analysis
            sensitivity.PerformAnalysis();
            
            // Save analysis to file
            string fileName = "primal_simplex_sensitivity.txt";
            string analysisText = sensitivity.GenerateVisualization();
            File.WriteAllText(fileName, analysisText);
            Console.WriteLine($"\nSensitivity analysis saved to: {fileName}");
        }
        
        private static void RunRevisedPrimalSimplexSensitivity(ParsedLinearProgrammingModel model)
        {
            Console.WriteLine("\n=== Revised Primal Simplex Sensitivity Analysis ===");
            var canonicalModel = model.ToCanonicalForm();
            var solver = new RevisedPrimalSimplexSolver();
            
            Console.WriteLine("Solving with Revised Primal Simplex...");
            var solution = solver.Solve(canonicalModel);
            
            if (!string.Equals(solution.Status, "Optimal", StringComparison.OrdinalIgnoreCase))
            {
                Console.WriteLine("Cannot perform sensitivity analysis - the problem does not have an optimal solution.");
                return;
            }
            
            // Get constraint types and RHS from the canonical model
            var constraintTypes = canonicalModel.ConstraintTypes;
            var rhsCoefficients = canonicalModel.RHSVector;
            
            // Initialize sensitivity analysis with required parameters
            var sensitivity = new SimplexSensitivityAnalysis(
                finalTableau: solution.OptimalTable,
                basisIndices: solution.BasisIndices,
                variableNames: solution.VariableNames,
                objectiveCoefficients: canonicalModel.ObjectiveCoefficients,
                rhsCoefficients: rhsCoefficients,
                constraintTypes: constraintTypes
            );
            
            // Perform sensitivity analysis
            sensitivity.PerformAnalysis();
            
            // Save analysis to file
            string fileName = "revised_primal_simplex_sensitivity.txt";
            string analysisText = sensitivity.GenerateVisualization();
            File.WriteAllText(fileName, analysisText);
            Console.WriteLine($"\nSensitivity analysis saved to: {fileName}");
        }
        
        private static void RunBranchAndBoundSensitivity(ParsedLinearProgrammingModel model)
        {
            Console.WriteLine("\n=== Branch and Bound Sensitivity Analysis ===");
            
            try
            {
                // First, we need to run the Branch and Bound solver to get a solution
                Console.WriteLine("Running Branch and Bound solver to get a solution...");
                var canonicalModel = model.ToCanonicalForm();
                var solver = new BranchAndBoundSolver();
                var solution = solver.Solve(model); // Pass the original model, not the canonical one
                
                if (solution == null)
                {
                    Console.WriteLine("Failed to find a solution using Branch and Bound.");
                    return;
                }
                
                // Now create the sensitivity analysis with the solution and original model
                var sensitivity = new BranchAndBoundSensitivityAnalysis(solution, model);
                
                // Perform the sensitivity analysis
                sensitivity.PerformAnalysis();
                
                // Generate and display the analysis results
                string analysisResults = sensitivity.GenerateVisualization();
                Console.WriteLine("\n=== Sensitivity Analysis Results ===");
                Console.WriteLine(analysisResults);
                
                // Save analysis to file
                string fileName = "branch_bound_sensitivity.txt";
                File.WriteAllText(fileName, analysisResults);
                Console.WriteLine($"\nSensitivity analysis saved to: {fileName}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nError during Branch and Bound sensitivity analysis: {ex.Message}");
                if (ex.InnerException != null)
                {
                    Console.WriteLine($"Inner exception: {ex.InnerException.Message}");
                }
            }
        }
        
        private static void RunCuttingPlaneSensitivity(ParsedLinearProgrammingModel model)
        {
            Console.WriteLine("\n=== Cutting Plane Sensitivity Analysis ===");
            
            try
            {
                // Convert the parsed model to canonical form first
                var parsedCanonicalModel = model.ToCanonicalForm();
                
                // Create a new CanonicalLinearProgrammingModel from the parsed canonical model
                var canonicalModel = new CanonicalLinearProgrammingModel
                {
                    CoefficientMatrix = parsedCanonicalModel.CoefficientMatrix,
                    RightHandSide = parsedCanonicalModel.RightHandSide,
                    ObjectiveCoefficients = parsedCanonicalModel.ObjectiveCoefficients,
                    ConstraintTypes = parsedCanonicalModel.ConstraintTypes,
                    OptimizationType = parsedCanonicalModel.OptimizationType,
                    VariableTypes = parsedCanonicalModel.VariableTypes
                };
                
                // Run the Cutting Plane solver with the canonical model
                Console.WriteLine("Running Cutting Plane solver to get a solution...");
                var solver = new CuttingPlane();
                var solution = solver.CuttingPlaneSolve(canonicalModel);
                
                if (solution == null)
                {
                    Console.WriteLine("Failed to find a solution using Cutting Plane method.");
                    return;
                }
                
                // Create the sensitivity analysis with the solution
                // Note: This is a simplified version - in a real implementation, you would need to provide
                // the LP relaxation solution and cutting plane history as well
                var sensitivity = new CuttingPlaneSensitivityAnalysis(
                    lpRelaxationSolution: solution, // In a real implementation, this should be the LP relaxation solution
                    finalSolution: solution,
                    originalModel: canonicalModel,
                    cuttingPlaneHistory: new List<double[,]>() // Empty history for now
                );
                
                // Perform the sensitivity analysis
                sensitivity.PerformAnalysis();
                
                // Generate and display the analysis results
                string analysisResults = sensitivity.GenerateVisualization();
                Console.WriteLine("\n=== Sensitivity Analysis Results ===");
                Console.WriteLine(analysisResults);
                
                // Save analysis to file
                string fileName = "cutting_plane_sensitivity.txt";
                File.WriteAllText(fileName, analysisResults);
                Console.WriteLine($"\nSensitivity analysis saved to {fileName}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nError during Cutting Plane sensitivity analysis: {ex.Message}");
                if (ex.InnerException != null)
                {
                    Console.WriteLine($"Inner exception: {ex.InnerException.Message}");
                }
            }
        }
        
        private static void PerformKnapsackSensitivityAnalysis(double[] weights, double[] values, double capacity, double maxValue, bool[] includedItems, string solutionPath)
        {
            try
            {
                Console.WriteLine("\n=== Knapsack Sensitivity Analysis ===");
                
                // Create a new knapsack solver instance
                var knapsack = new KnapsackSolver(weights, values, capacity);
                
                // Perform sensitivity analysis
                knapsack.PerformSensitivityAnalysis(includedItems, maxValue);
                
                // Save the sensitivity analysis to a file
                string analysisPath = Path.Combine(
                    Path.GetDirectoryName(solutionPath),
                    Path.GetFileNameWithoutExtension(solutionPath) + "_sensitivity.txt");
                
                File.WriteAllText(analysisPath, string.Join(Environment.NewLine, knapsack.IterationLog));
                Console.WriteLine($"\nSensitivity analysis saved to: {Path.GetFullPath(analysisPath)}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nError performing sensitivity analysis: {ex.Message}");
                Console.WriteLine(ex.StackTrace);
            }
        }
        
        private static void RunKnapsackSensitivity()
        {
            Console.WriteLine("\n=== Knapsack Problem Sensitivity Analysis ===");
            
            // Get the same input as the regular knapsack solver
            Console.WriteLine("First, provide the knapsack problem details:");
            
            int n;
            while (true)
            {
                Console.Write("Enter number of items: ");
                if (int.TryParse(Console.ReadLine(), out n) && n > 0)
                    break;
                Console.WriteLine("Please enter a positive integer.");
            }

            double[] weights = new double[n];
            double[] values = new double[n];

            for (int i = 0; i < n; i++)
            {
                Console.WriteLine($"\nItem {i + 1}:");
                
                while (true)
                {
                    Console.Write("  Weight: ");
                    if (double.TryParse(Console.ReadLine(), out double weight) && weight > 0)
                    {
                        weights[i] = weight;
                        break;
                    }
                    Console.WriteLine("  Please enter a positive number for weight.");
                }

                while (true)
                {
                    Console.Write("  Value: ");
                    if (double.TryParse(Console.ReadLine(), out double value) && value >= 0)
                    {
                        values[i] = value;
                        break;
                    }
                    Console.WriteLine("  Please enter a non-negative number for value.");
                }
            }

            double capacity;
            while (true)
            {
                Console.Write("\nEnter capacity of knapsack: ");
                if (double.TryParse(Console.ReadLine(), out capacity) && capacity > 0)
                    break;
                Console.WriteLine("Please enter a positive number for capacity.");
            }
            
            // Solve the knapsack problem first
            var solver = new KnapsackSolver(weights, values, capacity);
            var (maxValue, includedItems) = solver.Solve();
            
            // Perform and display sensitivity analysis
            var sensitivity = new KnapsackSensitivityAnalysis(weights, values, capacity, includedItems, maxValue);
            sensitivity.PerformAnalysis();
            
            // Save analysis to file
            string fileName = "knapsack_sensitivity.txt";
            File.WriteAllText(fileName, sensitivity.ToString());
            Console.WriteLine($"\nSensitivity analysis saved to {fileName}");
        }
    }

    /// <summary>
    /// Handles the generation of output files for different solvers
    /// </summary>
    public static class OutputFileGenerator
    {
        public static void GenerateSimplexOutput(
            ParsedLinearProgrammingModel model,
            LinearProgramming.Algorithms.PrimalSimplex.LinearProgramSolution solution,
            string algorithmName,
            string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                writer.WriteLine($"=== {algorithmName} Solution ===");
                writer.WriteLine($"Status: {solution.Status}");
                writer.WriteLine($"Objective Value: {solution.ObjectiveValue:F6}");
                
                if (solution.SolutionVector != null)
                {
                    writer.WriteLine("\nVariable Values:");
                    for (int i = 0; i < solution.SolutionVector.Length; i++)
                    {
                        string varName = i < model.Variables.Count ? model.Variables[i].Name : $"s{i - model.Variables.Count + 1}";
                        writer.WriteLine($"{varName} = {solution.SolutionVector[i]:F6}");
                    }
                }
                
                if (solution is LinearProgramming.Algorithms.PrimalSimplex.LinearProgramSolution simplexSolution && 
                    simplexSolution.BasisIndices != null)
                {
                    writer.WriteLine("\nBasis Variables:");
                    for (int i = 0; i < simplexSolution.BasisIndices.Length; i++)
                    {
                        int varIndex = simplexSolution.BasisIndices[i];
                        if (varIndex < model.Variables.Count)
                        {
                            string varName = model.Variables[varIndex].Name;
                            double varValue = simplexSolution.SolutionVector[varIndex];
                            writer.WriteLine($"{varName} = {varValue:F6}");
                        }
                    }
                }
            }
        }
        
        public static void GenerateIntegerSolutionOutput(
            ParsedLinearProgrammingModel model,
            dynamic solution,
            string algorithmName,
            string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                writer.WriteLine($"=== {algorithmName} Solution ===");
                writer.WriteLine($"Status: {solution.Status}");
                writer.WriteLine($"Objective Value: {solution.ObjectiveValue:F0}");
                
                if (solution.SolutionVector != null)
                {
                    writer.WriteLine("\nVariable Values:");
                    for (int i = 0; i < solution.SolutionVector.Length; i++)
                    {
                        string varName = i < model.Variables.Count ? model.Variables[i].Name : $"s{i - model.Variables.Count + 1}";
                        writer.WriteLine($"{varName} = {solution.SolutionVector[i]:F0}");
                    }
                }
            }
        }
        
        public static void GenerateBranchAndBoundOutput(
            ParsedLinearProgrammingModel model,
            LinearProgramming.Algorithms.BranchAndBound.BranchAndBoundSolution solution,
            string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                writer.WriteLine("=== Branch and Bound Solution ===");
                writer.WriteLine($"Status: {solution.OverallStatus}");
                writer.WriteLine($"Objective Value: {solution.BestObjectiveValue:F0}");
                writer.WriteLine($"Nodes Explored: {solution.TotalNodesExplored}");
                writer.WriteLine($"Nodes Fathomed: {solution.TotalNodesFathomed}");
                
                if (solution.BestSolution != null)
                {
                    writer.WriteLine("\nVariable Values:");
                    for (int i = 0; i < solution.BestSolution.Length; i++)
                    {
                        string varName = i < model.Variables.Count ? model.Variables[i].Name : $"s{i - model.Variables.Count + 1}";
                        writer.WriteLine($"{varName} = {solution.BestSolution[i]:F0}");
                    }
                }
            }
        }
        
        public static void GenerateKnapsackOutput(
            double[] weights,
            double[] values,
            double capacity,
            (double maxValue, bool[] includedItems) solution,
            string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                writer.WriteLine("=== Knapsack Solution ===");
                writer.WriteLine($"Capacity: {capacity:F2}");
                double totalWeight = 0;
                writer.WriteLine($"Maximum Value: {solution.maxValue:F2}");
                
                writer.WriteLine("\nSelected Items:");
                for (int i = 0; i < solution.includedItems.Length; i++)
                {
                    if (solution.includedItems[i])
                    {
                        totalWeight += weights[i];
                        writer.WriteLine($"  Item {i + 1}: Value = {values[i]:F2}, Weight = {weights[i]:F2}");
                    }
                }
                writer.WriteLine($"\nTotal Weight: {totalWeight:F2}");
                
                // List all items with selection status
                writer.WriteLine("\nAll Items:");
                for (int i = 0; i < weights.Length; i++)
                {
                    string selected = solution.includedItems[i] ? "[X]" : "[ ]";
                    writer.WriteLine($"{selected} Item {i + 1}: Value = {values[i]:F2}, Weight = {weights[i]:F2}");
                }
            }
        }
    }
}
