using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.Knapsack
{
    public class KnapsackNode
    {
        public int Level { get; set; }
        public double Value { get; set; }
        public double Weight { get; set; }
        public double Bound { get; set; }
        public bool[] Included { get; set; }
        public string Path { get; set; }

        public KnapsackNode(int n)
        {
            Included = new bool[n];
            Path = "";
        }
    }

    public class KnapsackSolver
    {
        private double[] weights;
        private double[] values;
        private double capacity;
        private int n;
        private KnapsackNode bestNode;
        private List<string> iterationLog;
        private const double Epsilon = 1e-10;

        private static bool AreEqual(double a, double b, double epsilon = 1e-10)
        {
            return Math.Abs(a - b) < epsilon;
        }

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
        
        public KnapsackSolver(double[] weights, double[] values, double capacity)
        {
            this.weights = weights ?? throw new ArgumentNullException(nameof(weights));
            this.values = values ?? throw new ArgumentNullException(nameof(values));
            this.capacity = capacity;
            this.n = weights.Length;
            this.iterationLog = new List<string>();
            
            // Validate inputs and throw detailed exception if invalid
            string validationError = ValidateInputs();
            if (validationError != null)
            {
                throw new ArgumentException($"Invalid knapsack problem: {validationError}");
            }
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
        /// Prepares the knapsack problem by sorting items and checking for trivial solutions
        /// </summary>
        private bool PrepareProblem(out double maxValue, out bool[] includedItems)
        {
            maxValue = 0;
            includedItems = new bool[n];
            
            // Check for trivial case: no items or zero capacity
            if (n == 0 || capacity <= 0)
            {
                return true;
            }
            
            // Check for trivial case: only one item that fits
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
                for (int i = 0; i < n; i++) includedItems[i] = true;
                return true;
            }
            
            return false;
        }
        
        public (double maxValue, bool[] includedItems) Solve()
        {
            // Check for trivial solutions first
            if (PrepareProblem(out double trivialValue, out bool[] trivialSolution))
            {
                LogTrivialSolution(trivialValue, trivialSolution);
                return (trivialValue, trivialSolution);
            }
            
            // Sort items by value/weight ratio in descending order
            var items = Enumerable.Range(0, n)
                .OrderByDescending(i => values[i] / weights[i])
                .ToList();

            // Rearrange weights and values based on sorted order
            weights = items.Select(i => weights[i]).ToArray();
            values = items.Select(i => values[i]).ToArray();

            var queue = new Queue<KnapsackNode>();
            var root = new KnapsackNode(n) { Level = -1, Value = 0, Weight = 0 };
            root.Bound = CalculateBound(root);
            queue.Enqueue(root);

            while (queue.Count > 0)
            {
                var currentNode = queue.Dequeue();
                
                // If there's no chance of doing better, skip this node
                if (bestNode != null && currentNode.Bound <= bestNode.Value)
                {
                    LogNode(currentNode, "Fathomed (bound <= best value)");
                    continue;
                }

                // If we've reached the end, check if this is the best solution
                if (currentNode.Level == n - 1)
                {
                    if (currentNode.Value > (bestNode?.Value ?? 0) && currentNode.Weight <= capacity)
                    {
                        bestNode = currentNode;
                        LogNode(bestNode, "New best solution found!");
                    }
                    continue;
                }

                // Create child node where we include the next item
                var includeNode = new KnapsackNode(n)
                {
                    Level = currentNode.Level + 1,
                    Weight = currentNode.Weight + weights[currentNode.Level + 1],
                    Path = currentNode.Path + "1"
                };
                Array.Copy(currentNode.Included, includeNode.Included, n);
                includeNode.Included[includeNode.Level] = true;
                includeNode.Value = currentNode.Value + values[includeNode.Level];
                includeNode.Bound = CalculateBound(includeNode);

                // Create child node where we exclude the next item
                var excludeNode = new KnapsackNode(n)
                {
                    Level = currentNode.Level + 1,
                    Weight = currentNode.Weight,
                    Value = currentNode.Value,
                    Path = currentNode.Path + "0"
                };
                Array.Copy(currentNode.Included, excludeNode.Included, n);
                excludeNode.Included[excludeNode.Level] = false;
                excludeNode.Bound = CalculateBound(excludeNode);

                // Add nodes to queue if they're promising
                if (includeNode.Weight <= capacity)
                {
                    queue.Enqueue(includeNode);
                    LogNode(includeNode, "Include item " + (includeNode.Level + 1));
                }
                else
                {
                    LogNode(includeNode, "Fathomed (over capacity)");
                }

                if (excludeNode.Bound > (bestNode?.Value ?? 0))
                {
                    queue.Enqueue(excludeNode);
                    LogNode(excludeNode, "Exclude item " + (excludeNode.Level + 1));
                }
                else
                {
                    LogNode(excludeNode, "Fathomed (bound <= best value)");
                }
            }

            // Map back to original item order if we sorted
            var result = new bool[n];
            if (bestNode != null)
            {
                for (int i = 0; i < n; i++)
                {
                    result[i] = bestNode.Included[i];
                }
                return (bestNode.Value, result);
            }

            return (0, new bool[n]);
        }

        private void LogNode(KnapsackNode node, string action)
        {
            var log = new StringBuilder();
            log.AppendLine(OutputFormatter.CreateHeader($"Node: {node.Path} - {action}"));
            log.AppendLine(OutputFormatter.FormatKeyValue("Level", node.Level + 1));
            log.AppendLine(OutputFormatter.FormatKeyValue("Current Value", node.Value));
            log.AppendLine(OutputFormatter.FormatKeyValue("Current Weight", $"{node.Weight} / {capacity}"));
            log.AppendLine(OutputFormatter.FormatKeyValue("Upper Bound", node.Bound));
            
            if (bestNode != null)
            {
                log.AppendLine(OutputFormatter.FormatKeyValue("Best Value So Far", bestNode.Value));
            }
            
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