using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.Knapsack
{
    /// <summary>
    /// Represents a node in the branch and bound search tree for the knapsack problem
    /// </summary>
    public class KnapsackNode
    {
        public int Level { get; set; }
        public double Value { get; set; }
        public double Weight { get; set; }
        public double Bound { get; set; }
        public bool[] Included { get; set; }
        public string Path { get; set; }
        public int? ParentId { get; set; }

        public KnapsackNode(int n)
        {
            Included = new bool[n];
            Path = "";
            ParentId = null;
        }
        
        public KnapsackNode Clone()
        {
            return new KnapsackNode(Included.Length)
            {
                Level = this.Level,
                Value = this.Value,
                Weight = this.Weight,
                Bound = this.Bound,
                Included = (bool[])this.Included.Clone(),
                Path = this.Path,
                ParentId = this.ParentId
            };
        }
    }

    /// <summary>
    /// Solves the 0/1 knapsack problem using branch and bound algorithm
    /// </summary>
    public class KnapsackSolver
    {
        private double[] weights;
        private double[] values;
        private readonly double capacity;
        private const double Epsilon = 1e-10;
        private int n;
        private KnapsackNode bestNode;
        private readonly List<string> iterationLog = new();
        private readonly List<KnapsackNode> searchTree = new();

        public KnapsackSolver(double[] weights, double[] values, double capacity)
        {
            this.weights = weights ?? throw new ArgumentNullException(nameof(weights));
            this.values = values ?? throw new ArgumentNullException(nameof(values));
            this.capacity = capacity > 0 ? capacity : throw new ArgumentException("Capacity must be positive", nameof(capacity));
            
            if (weights.Length != values.Length)
            {
                throw new ArgumentException("Number of weights must match number of values");
            }
            
            this.n = weights.Length;
            this.iterationLog = new List<string>();
            this.searchTree = new List<KnapsackNode>();
        }
        
        /// <summary>
        /// Gets the iteration log containing detailed information about the solving process
        /// </summary>
        public IReadOnlyList<string> IterationLog => iterationLog.AsReadOnly();
        
        /// <summary>
        /// Gets the search tree generated during the solving process
        /// </summary>
        public IReadOnlyList<KnapsackNode> SearchTree => searchTree.AsReadOnly();
        
        /// <summary>
        /// Gets the best node found during the search
        /// </summary>
        public KnapsackNode BestNode => bestNode;

        private static bool AreEqual(double a, double b, double epsilon = 1e-10)
        {
            return Math.Abs(a - b) < epsilon;
        }

        /// <summary>
        /// Validates the knapsack problem inputs and checks for trivially infeasible cases
        /// </summary>
        /// <summary>
        /// Validates the knapsack problem inputs and checks for trivially infeasible cases
        /// </summary>
        private string ValidateInputs()
        {
            var issues = new List<string>();
            
            // Check for null or empty inputs
            if (weights == null || values == null)
                return "Weights and values arrays cannot be null";
                
            if (weights.Length != values.Length)
                return "Length of weights and values arrays must be equal";

            if (weights.Length == 0)
                return "At least one item must be provided";
                
            // Check for negative or invalid values and potential overflow
            double totalWeight = 0;
            double totalValue = 0;
            
            for (int i = 0; i < n; i++)
            {
                // Check for invalid numbers
                if (double.IsNaN(weights[i]) || double.IsInfinity(weights[i]))
                    issues.Add($"Item {i + 1} has invalid weight: {weights[i]}");
                    
                if (double.IsNaN(values[i]) || double.IsInfinity(values[i]))
                    issues.Add($"Item {i + 1} has invalid value: {values[i]}");
                
                // Check for negative values
                if (weights[i] < 0)
                    issues.Add($"Item {i + 1} has negative weight ({weights[i]})");
                    
                if (values[i] < 0)
                    issues.Add($"Item {i + 1} has negative value ({values[i]})");
                
                // Check for potential overflow in totals
                try
                {
                    checked
                    {
                        totalWeight += weights[i];
                        totalValue += values[i];
                    }
                }
                catch (OverflowException)
                {
                    issues.Add("Calculation would result in numeric overflow. Please use smaller input values.");
                    break;
                }
                
                // Check for zero or extremely small weights to avoid division by zero
                if (weights[i].IsZero(Epsilon * 10))
                {
                    issues.Add($"Item {i + 1} has zero or negligible weight. This may cause division by zero.");
                }
            }
            
            // Check capacity
            if (capacity < 0)
                issues.Add("Knapsack capacity cannot be negative");
            else if (capacity == 0)
                issues.Add("Knapsack capacity is zero - no items can be included");
                
            // Check if any items can fit
            if (capacity > 0 && weights.All(w => w > capacity))
                issues.Add("No items can fit in the knapsack (all items exceed capacity)");
                
            return issues.Count > 0 ? string.Join("\n", issues) : null;
        }
        
        /// <summary>
        /// Calculates the upper bound for a node using the fractional knapsack approach
        /// </summary>
        private double CalculateBound(KnapsackNode node)
        {
            try
            {
                if (node.Weight > capacity + Epsilon)
                {
                    LogNode(node, "Fathomed (knapsack full)");
                    return 0;
                }

                double bound = node.Value;
                int j = node.Level + 1;
                double totalWeight = node.Weight;

                // Greedily add items until we can't fit more
                while (j < n && (totalWeight + weights[j] <= capacity + Epsilon))
                {
                    totalWeight = (totalWeight + weights[j]).RoundToEpsilon(Epsilon);
                    bound = (bound + values[j]).RoundToEpsilon(Epsilon);
                    j++;
                }

                // Add fraction of next item if there's still space
                if (j < n && !weights[j].IsZero(Epsilon))
                {
                    double remainingCapacity = (capacity - totalWeight).Clamp(0, capacity);
                    double fraction = remainingCapacity / weights[j];
                    bound += fraction * values[j];
                    bound = bound.RoundToEpsilon(Epsilon);
                }

                return bound;
            }
            catch (Exception ex) when (ex is OverflowException || ex is DivideByZeroException)
            {
                // Log the error and return a safe bound
                LogNode(node, $"Error calculating bound: {ex.Message}");
                return node.Value; // Return at least the current value
            }
        }

        /// <summary>
        /// Logs information about a node being processed
        /// </summary>
        private void LogNode(KnapsackNode node, string action)
        {
            var log = new StringBuilder();
            log.AppendLine(OutputFormatter.CreateHeader($"Node: {node.Path} - {action}"));
            log.AppendLine(OutputFormatter.FormatKeyValue("Level", node.Level + 1));
            log.AppendLine(OutputFormatter.FormatKeyValue("Current Value", node.Value));
            log.AppendLine(OutputFormatter.FormatKeyValue("Current Weight", $"{node.Weight} / {capacity}"));
            log.AppendLine(OutputFormatter.FormatKeyValue("Upper Bound", node.Bound));
            
            var included = Enumerable.Range(0, n)
                .Where(i => node.Included[i])
                .Select(i => $"Item {i + 1}");
                
            log.AppendLine(OutputFormatter.FormatKeyValue("Included Items", 
                included.Any() ? string.Join(", ", included) : "None"));
                
            iterationLog.Add(log.ToString());
        }
        
        /// <summary>
        /// Logs the current best solution
        /// </summary>
        private void LogCurrentBest()
        {
            if (bestNode == null) return;
            
            var log = new StringBuilder();
            log.AppendLine(OutputFormatter.CreateHeader("CURRENT BEST SOLUTION"));
            log.AppendLine(OutputFormatter.FormatKeyValue("Value", bestNode.Value));
            log.AppendLine(OutputFormatter.FormatKeyValue("Weight", $"{bestNode.Weight} / {capacity}"));
            
            var included = Enumerable.Range(0, n)
                .Where(i => bestNode.Included[i])
                .Select(i => $"Item {i + 1}");
                
            log.AppendLine(OutputFormatter.FormatKeyValue("Items", 
                included.Any() ? string.Join(", ", included) : "None"));
                
            iterationLog.Add(log.ToString());
        }

        /// <summary>
        /// Solves the 0/1 knapsack problem using branch and bound algorithm
        /// </summary>
        /// <param name="enableLogging">Whether to enable detailed logging of the solving process</param>
        /// <returns>A tuple containing the maximum value and a boolean array indicating included items</returns>
        public (double maxValue, bool[] includedItems) Solve(bool enableLogging = true)
        {
            // Reset state
            iterationLog.Clear();
            searchTree.Clear();
            bestNode = null;
            
            // Check for trivial solutions first
            if (PrepareProblem(out double trivialValue, out bool[] trivialSolution))
            {
                if (enableLogging)
                {
                    LogTrivialSolution(trivialValue, trivialSolution);
                }
                return (trivialValue, trivialSolution);
            }
            
            // Sort items by value/weight ratio in descending order
            var items = Enumerable.Range(0, n)
                .OrderByDescending(i => values[i] / (weights[i] + Epsilon))
                .ToArray();

            // Reorder arrays based on sorted items
            this.weights = items.Select(i => weights[i]).ToArray();
            this.values = items.Select(i => values[i]).ToArray();

            // Initialize best solution
            bestNode = new KnapsackNode(n)
            {
                Level = -1,
                Value = 0,
                Weight = 0,
                Bound = 0,
                Included = new bool[n]
            };

            // Create root node
            var root = new KnapsackNode(n)
            {
                Level = -1,
                Value = 0,
                Weight = 0,
                Bound = CalculateBound(new KnapsackNode(n)),
                Included = new bool[n],
                Path = "Root"
            };
            
            // Add root to search tree
            searchTree.Add(root);
            
            // Use a stack for DFS
            var nodeStack = new Stack<KnapsackNode>();
            nodeStack.Push(root);

            while (nodeStack.Count > 0)
            {
                var currentNode = nodeStack.Pop();
                
                // If bound is worse than best, skip this node
                if (currentNode.Bound <= bestNode.Value + Epsilon)
                {
                    if (enableLogging)
                    {
                        LogNode(currentNode, "Fathomed (bound <= best)");
                    }
                    continue;
                }

                // If we've reached the end or found a better solution, update best solution
                if (currentNode.Level == n - 1 || 
                    (currentNode.Value > bestNode.Value + Epsilon && currentNode.Weight <= capacity + Epsilon))
                {
                    if (currentNode.Value > bestNode.Value + Epsilon)
                    {
                        bestNode = currentNode.Clone();
                        if (enableLogging)
                        {
                            LogNode(bestNode, "New best solution found");
                            LogCurrentBest();
                        }
                    }
                    
                    if (currentNode.Level == n - 1)
                    {
                        continue;
                    }
                }

                int nextLevel = currentNode.Level + 1;

                // Right child (exclude next item)
                var rightNode = new KnapsackNode(n)
                {
                    Level = nextLevel,
                    Value = currentNode.Value,
                    Weight = currentNode.Weight,
                    Included = (bool[])currentNode.Included.Clone(),
                    Path = $"{currentNode.Path} -> Exclude {nextLevel + 1}",
                    ParentId = searchTree.IndexOf(currentNode)
                };
                rightNode.Bound = CalculateBound(rightNode);
                
                // Left child (include next item)
                var leftNode = new KnapsackNode(n)
                {
                    Level = nextLevel,
                    Value = currentNode.Value + values[nextLevel],
                    Weight = currentNode.Weight + weights[nextLevel],
                    Included = (bool[])currentNode.Included.Clone(),
                    Path = $"{currentNode.Path} -> Include {nextLevel + 1}",
                    ParentId = searchTree.IndexOf(currentNode)
                };
                leftNode.Included[nextLevel] = true;
                leftNode.Bound = CalculateBound(leftNode);

                // Add right child first (LIFO for DFS)
                if (rightNode.Weight <= capacity + Epsilon)
                {
                    rightNode.ParentId = searchTree.IndexOf(currentNode);
                    searchTree.Add(rightNode);
                    nodeStack.Push(rightNode);
                    if (enableLogging)
                    {
                        LogNode(rightNode, "Created right child (exclude)");
                    }
                }

                // Add left child if it's feasible
                if (leftNode.Weight <= capacity + Epsilon)
                {
                    leftNode.ParentId = searchTree.IndexOf(currentNode);
                    searchTree.Add(leftNode);
                    nodeStack.Push(leftNode);
                    if (enableLogging)
                    {
                        LogNode(leftNode, "Created left child (include)");
                    }
                }
                
                // Log current best solution
                if (enableLogging)
                {
                    LogCurrentBest();
                }
            }

            // Reconstruct the solution in original item order
            var solution = new bool[n];
            for (int i = 0; i < n; i++)
            {
                if (bestNode.Included[i])
                {
                    solution[items[i]] = true;
                }
            }

            return (bestNode.Value, solution);
        }

        /// <summary>
        /// Logs information about a trivial solution
        /// </summary>
        private void LogTrivialSolution(double value, bool[] solution)
        {
            var log = new StringBuilder();
            log.AppendLine(OutputFormatter.CreateHeader("TRIVIAL SOLUTION FOUND"));
            
            if (n == 0)
            {
                log.AppendLine("No items to consider");
            }
            else if (capacity <= 0)
            {
                log.AppendLine("Knapsack capacity is zero");
            }
            else if (n == 1)
            {
                log.AppendLine("Only one item to consider");
                log.AppendLine(OutputFormatter.FormatKeyValue("Item fits", solution[0] ? "Yes" : "No"));
            }
            else
            {
                log.AppendLine("All items fit in the knapsack");
            }
            
            log.AppendLine(OutputFormatter.FormatKeyValue("Total Value", value));
            log.AppendLine(OutputFormatter.FormatKeyValue("Total Weight", 
                solution.Zip(weights, (inc, w) => inc ? w : 0).Sum()));
                
            iterationLog.Add(log.ToString());
        }

        /// <summary>
        /// Prepares the knapsack problem by sorting items and checking for trivial solutions
        /// </summary>
        private bool PrepareProblem(out double maxValue, out bool[] includedItems)
        {
            maxValue = 0;
            includedItems = new bool[n];
            
            // Check for empty problem
            if (n == 0)
            {
                return true;
            }
            
            // Check for zero capacity
            if (capacity <= 0)
            {
                return true;
            }
            
            // Check for single item
            if (n == 1)
            {
                if (weights[0] <= capacity)
                {
                    maxValue = values[0];
                    includedItems[0] = true;
                }
                return true;
            }
            
            // Check if all items can fit
            double totalWeight = weights.Sum();
            if (totalWeight <= capacity)
            {
                maxValue = values.Sum();
                for (int i = 0; i < n; i++) 
                    includedItems[i] = true;
                return true;
            }

            return false; // No trivial solution found
        }

        // Create root node in the Solve method
        private KnapsackNode CreateRootNode()
        {
            return new KnapsackNode(n)
            {
                Level = -1,
                Value = 0,
                Weight = 0,
                Bound = CalculateBound(new KnapsackNode(n)),
                Included = new bool[n],
                Path = "Root"
            };
        }

        private void LogItems(StringBuilder log, KnapsackNode node, double[] values, double[] weights)
        {
            log.Append("Items: ");
            bool anyIncluded = false;
            for (int i = 0; i <= node.Level; i++)
            {
                if (node.Included[i])
                {
                    log.Append($"{i + 1} " + $"(w:{weights[i]},v:{values[i]}) ");
                    anyIncluded = true;
                }
            }
            
            if (!anyIncluded)
            {
                log.Append("None");
            }
            iterationLog.Add(log.ToString());
        }

        /// <summary>
        /// Performs sensitivity analysis on the knapsack solution
        /// </summary>
        /// <param name="solution">The solution array where true indicates the item is included</param>
        /// <param name="optimalValue">The optimal value of the knapsack solution</param>
        public void PerformSensitivityAnalysis(bool[] solution, double optimalValue)
        {
            try
            {
                Console.WriteLine("\n" + new string('=', 80));
                Console.WriteLine("KNAPSACK SENSITIVITY ANALYSIS");
                Console.WriteLine(new string('=', 80));

                var sensitivity = new KnapsackSensitivityAnalysis(
                    weights, 
                    values, 
                    capacity, 
                    solution, 
                    optimalValue);
                    
                sensitivity.PerformAnalysis();
                
                Console.WriteLine("\n" + new string('=', 80));
                Console.WriteLine("END OF SENSITIVITY ANALYSIS");
                Console.WriteLine(new string('=', 80));
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\n⚠ Error during sensitivity analysis: {ex.Message}");
                Console.WriteLine(ex.StackTrace);
            }
        }

        /// <summary>
        /// Performs sensitivity analysis on the best found solution
        /// </summary>
        public void PerformSensitivityAnalysis()
        {
            if (bestNode == null)
            {
                Console.WriteLine("No solution available for sensitivity analysis.");
                return;
            }
            PerformSensitivityAnalysis(bestNode.Included, bestNode.Value);
        }

        public void PrintIterations()
        {
            // Display the iteration log if available
            if (iterationLog.Count > 0)
            {
                Console.WriteLine(OutputFormatter.CreateHeader("KNAPSACK - BRANCH AND BOUND ITERATIONS"));
                foreach (var log in iterationLog)
                {
                    Console.WriteLine(log);
                }
            }

            // Display the best solution found
            if (bestNode != null)
            {
                Console.WriteLine(OutputFormatter.CreateHeader("BEST SOLUTION FOUND"));
                Console.WriteLine(OutputFormatter.FormatKeyValue("Status", "Optimal"));
                Console.WriteLine(OutputFormatter.FormatKeyValue("Total Value", bestNode.Value));
                Console.WriteLine(OutputFormatter.FormatKeyValue("Total Weight", $"{bestNode.Weight} / {capacity}"));
                
                // Calculate and display capacity utilization
                double utilization = (bestNode.Weight / capacity) * 100;
                Console.WriteLine(OutputFormatter.FormatKeyValue("Capacity Used", $"{utilization:F2}%"));
                
                // Display included items with details
                Console.WriteLine("\nINCLUDED ITEMS:");
                Console.WriteLine(string.Format("{0,-6} {1,-10} {2,-10} {3}", "Item", "Value", "Weight", "Value/Weight"));
                Console.WriteLine(new string('-', 45));
                
                for (int i = 0; i < n; i++)
                {
                    if (bestNode.Included[i])
                    {
                        double vwRatio = values[i] / weights[i];
                        Console.WriteLine($"{i + 1,-6} {values[i],-10:F2} {weights[i],-10:F2} {vwRatio,-10:F4}");
                    }
                }
                
                // Display excluded items for reference
                var excludedItems = Enumerable.Range(0, n)
                    .Where(i => !bestNode.Included[i])
.OrderByDescending(i => values[i] / weights[i])
                    .ToList();
                    
                if (excludedItems.Any())
                {
                    Console.WriteLine("\nEXCLUDED ITEMS (by value/weight ratio):");
                    Console.WriteLine(string.Format("{0,-6} {1,-10} {2,-10} {3}", "Item", "Value", "Weight", "Value/Weight"));
                    Console.WriteLine(new string('-', 45));
                    
                    foreach (var i in excludedItems.Take(10)) // Show top 10 excluded items
                    {
                        double vwRatio = values[i] / weights[i];
                        Console.WriteLine($"{i + 1,-6} {values[i],-10:F2} {weights[i],-10:F2} {vwRatio,-10:F4}");
                    }
                    
                    if (excludedItems.Count > 10)
                    {
                        Console.WriteLine($"... and {excludedItems.Count - 10} more items");
                    }
                }
                
                Console.WriteLine("\n" + new string('=', 50));
            }
        }
    }

    public static class Program
    {
        public static void Main(string[] args)
        {
            try
            {
                Console.WriteLine(OutputFormatter.CreateHeader("KNAPSACK PROBLEM - BRANCH AND BOUND SOLVER"));
                
                Console.WriteLine("Enter the number of items:");
                if (!int.TryParse(Console.ReadLine(), out int n) || n <= 0)
                {
                    throw new ArgumentException("Number of items must be a positive integer.");
                }

                double[] weights = new double[n];
                double[] values = new double[n];

                for (int i = 0; i < n; i++)
                {
                    Console.WriteLine($"\nItem {i + 1}:");
                    
                    Console.Write("  Weight: ");
                    if (!double.TryParse(Console.ReadLine(), out weights[i]) || weights[i] < 0)
                    {
                        throw new ArgumentException("Weight must be a non-negative number.");
                    }

                    Console.Write("  Value: ");
                    if (!double.TryParse(Console.ReadLine(), out values[i]) || values[i] < 0)
                    {
                        throw new ArgumentException("Value must be a non-negative number.");
                    }
                }

                Console.WriteLine("\nEnter the capacity of the knapsack:");
                if (!double.TryParse(Console.ReadLine(), out double capacity) || capacity <= 0)
                {
                    throw new ArgumentException("Capacity must be a positive number.");
                }

                var solver = new KnapsackSolver(weights, values, capacity);
                var (maxValue, includedItems) = solver.Solve();
                
                solver.PrintIterations();

                Console.WriteLine("\n" + OutputFormatter.CreateHeader("FINAL RESULT"));
                Console.WriteLine(OutputFormatter.FormatKeyValue("Maximum value", maxValue));
                Console.Write("Items included: ");
                for (int i = 0; i < includedItems.Length; i++)
                {
                    if (includedItems[i])
                    {
                        Console.Write($"{i + 1} ");
                    }
                }
                Console.WriteLine();
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nError: {ex.Message}");
            }
        }
    }
}