using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.Knapsack
{
    /// <summary>
    /// Provides sensitivity analysis for the Knapsack problem solution
    /// </summary>
    public class KnapsackSensitivityAnalysis : BaseSensitivityAnalysis, ISensitivityAnalysis
    {
        private readonly double[] _weights;
        private readonly double[] _values;
        private readonly double _capacity;
        private readonly bool[] _solution;
        private readonly double _optimalValue;
        private readonly int _itemCount;
        
        private void PrintSectionHeader(string title)
        {
            Console.WriteLine("\n" + new string('=', 80));
            Console.WriteLine(title);
            Console.WriteLine(new string('-', title.Length));
        }
        
        private void PrintSectionFooter()
        {
            Console.WriteLine(new string('=', 80) + "\n");
        }
        
        // Visualization data storage
        private Dictionary<string, (double min, double current, double max)> _objectiveRanges = new();
        private Dictionary<string, (double min, double current, double max)> _rhsRanges = new();
        private Dictionary<string, double> _shadowPrices = new();
        private Dictionary<string, double> _reducedCosts = new();
        // Removed _totalWeight field to avoid redundancy
        protected bool IsZero(double value)
        {
            return Math.Abs(value) < Epsilon;
        }
        

        /// <summary>
        /// Initializes a new instance of the KnapsackSensitivityAnalysis class
        /// </summary>
        public KnapsackSensitivityAnalysis(double[] weights, double[] values, double capacity, bool[] solution, double optimalValue)
            : base(
                Enumerable.Range(1, weights?.Length ?? 0).Select(i => $"Item {i}").ToArray(),
                new[] { ConstraintType.LessThanOrEqual },
                weights?.Length ?? 0,
                1)  // Single capacity constraint
        {
            _weights = weights ?? throw new ArgumentNullException(nameof(weights));
            _values = values ?? throw new ArgumentNullException(nameof(values));
            _capacity = capacity;
            _solution = solution ?? throw new ArgumentNullException(nameof(solution));
            _optimalValue = optimalValue;
            _itemCount = weights.Length;
            // Remove _totalWeight field since we'll calculate it on demand

            if (values.Length != _itemCount)
                throw new ArgumentException("Number of values must match number of weights");
            if (solution.Length != _itemCount)
                throw new ArgumentException("Solution array length must match number of items");
        }

        /// <summary>
        /// Performs complete sensitivity analysis on the Knapsack solution
        /// </summary>
        /// <summary>
        /// Performs complete sensitivity analysis on the Knapsack solution
        /// </summary>
        public override void PerformAnalysis()
        {
            PrintSectionHeader("KNAPSACK SENSITIVITY ANALYSIS");
            
            // Clear previous visualization data
            _objectiveRanges.Clear();
            _rhsRanges.Clear();
            _shadowPrices.Clear();
            _reducedCosts.Clear();

            // 1. Basic solution information
            DisplaySolutionInfo();

            // 2. Capacity sensitivity
            AnalyzeCapacitySensitivity();

            // 3. Item value sensitivity
            AnalyzeValueSensitivity();

            // 4. Weight sensitivity
            AnalyzeWeightSensitivity();

            // 5. Critical items analysis
            AnalyzeCriticalItems();
            
            // 6. Display visualizations
            DisplayVisualizations();
            
            PrintSectionFooter();
        }

        private void DisplaySolutionInfo()
        {
            PrintSubsectionHeader("SOLUTION INFORMATION");
            Console.WriteLine($"Optimal Value: {_optimalValue:F2}");
            double totalWeight = CalculateTotalWeight(_solution, _weights);
            Console.WriteLine($"Total Weight: {totalWeight:F2} / {_capacity:F2} ({totalWeight / _capacity * 100:F1}% of capacity)");
            Console.WriteLine($"Items in Knapsack: {_solution.Count(x => x)} of {_itemCount}");
            
            // Store RHS range for capacity constraint
            _rhsRanges["Capacity"] = (0, _capacity, double.PositiveInfinity);
            
            // Store objective coefficients and reduced costs
            for (int i = 0; i < _itemCount; i++)
            {
                string varName = $"Item {i + 1}";
                _objectiveRanges[varName] = (0, _values[i], _values[i]); // Will be updated in AnalyzeValueSensitivity
                _reducedCosts[varName] = _solution[i] ? 0 : _values[i]; // Simplified reduced cost
            }
            
            Console.WriteLine("\nSelected Items:");
            for (int i = 0; i < _itemCount; i++)
            {
                if (_solution[i])
                {
                    Console.WriteLine($"- Item {i + 1}: Value = {_values[i]:F2}, Weight = {_weights[i]:F2}, " +
                                    $"Value/Weight = {_values[i] / _weights[i]:F2}");
                }
            }
        }

        #region ISensitivityAnalysis Implementation
        
        /// <inheritdoc />
        public override (double lowerBound, double upperBound) GetRHSRange(int constraintIndex)
        {
            if (constraintIndex != 0)
                throw new ArgumentOutOfRangeException(nameof(constraintIndex), "Knapsack has only one constraint (capacity)");
                
            if (_rhsRanges.TryGetValue("Capacity", out var range))
                return (range.min, range.max);
            return (0, double.PositiveInfinity);
        }
        
        /// <inheritdoc />
        public override (double lowerBound, double upperBound) GetObjectiveCoefficientRange(int variableIndex)
        {
            if (variableIndex < 0 || variableIndex >= _itemCount)
                throw new ArgumentOutOfRangeException(nameof(variableIndex));
                
            string varName = $"Item {variableIndex + 1}";
            if (_objectiveRanges.TryGetValue(varName, out var range))
                return (range.min, range.max);
            return (0, double.PositiveInfinity);
        }
        
        /// <inheritdoc />
        public override double GetReducedCost(int variableIndex)
        {
            if (variableIndex < 0 || variableIndex >= _itemCount)
                throw new ArgumentOutOfRangeException(nameof(variableIndex));
                
            string varName = $"Item {variableIndex + 1}";
            return _reducedCosts.TryGetValue(varName, out var cost) ? cost : 0;
        }
        
        /// <inheritdoc />
        public override double GetShadowPrice(int constraintIndex)
        {
            if (constraintIndex != 0)
                throw new ArgumentOutOfRangeException(nameof(constraintIndex), "Knapsack has only one constraint (capacity)");
                
            return _shadowPrices.TryGetValue("Capacity", out var price) ? price : 0;
        }
        
        /// <inheritdoc />
        public override Dictionary<string, double> GetSolutionValues()
        {
            return _variableNames
                .Select((name, index) => new { name, value = _solution[index] ? 1.0 : 0.0 })
                .ToDictionary(x => x.name, x => x.value);
        }
        
        #endregion
        
        #region Visualization
        
        private void DisplayVisualizations()
        {
            PrintSubsectionHeader("VISUALIZATION");
            
            // Display objective coefficient ranges
            Console.WriteLine("\nObjective Coefficient Ranges (per item):");
            var objectiveData = _objectiveRanges
                .OrderByDescending(kv => kv.Value.current)
                .ToDictionary(kv => kv.Key, kv => kv.Value);
            Console.WriteLine(SensitivityVisualizer.CreateBarChart("Objective Coefficient Ranges", objectiveData, 60));
            
            // Display shadow prices (for capacity constraint)
            if (_shadowPrices.Any())
            {
                Console.WriteLine("\nShadow Price (Value per Unit of Capacity):");
                var shadowPriceData = new Dictionary<string, double>
                {
                    { "Capacity", _shadowPrices.Values.First() }
                };
                Console.WriteLine(SensitivityVisualizer.CreateHistogram(
                    "Shadow Price",
                    shadowPriceData));
            }
            
            // Display reduced costs
            Console.WriteLine("\nReduced Costs (Potential Value if Item is Added/Removed):");
            var reducedCosts = _reducedCosts
                .OrderByDescending(kv => Math.Abs(kv.Value))
                .Where(kv => !IsZero(kv.Value)) // Only show non-zero reduced costs
                .ToDictionary(kv => kv.Key, kv => kv.Value);
            
            if (reducedCosts.Any())
            {
                Console.WriteLine(SensitivityVisualizer.CreateHistogram(
                    "Reduced Costs",
                    reducedCosts));
            }
            else
            {
                Console.WriteLine("  No reduced costs (all items are at their optimal inclusion/exclusion)");
            }
        }
        
        #endregion
        
        private void AnalyzeCapacitySensitivity()
        {
            Console.WriteLine("\nCAPACITY SENSITIVITY");
            Console.WriteLine(new string('-', 40));

            // Calculate slack capacity using safe comparison
            double totalWeight = CalculateTotalWeight(_solution, _weights);
            double slack = _capacity - totalWeight;
            Console.WriteLine($"Current Slack: {slack:F6}");

            // Calculate minimum capacity reduction before solution becomes infeasible
            double minReduction = double.PositiveInfinity;
            for (int i = 0; i < _itemCount; i++)
            {
                if (_solution[i])
                {
                    minReduction = Math.Min(minReduction, _weights[i]);
                }
            }
            
            // Check if the knapsack is full (within numerical tolerance)
            bool isFull = Math.Abs(slack) < Epsilon * 10;
            Console.WriteLine($"Knapsack is {(isFull ? "full" : "not full")} (slack: {slack:E2})");
            
            if (isFull)
            {
                Console.WriteLine("Minimum capacity increase to allow any additional item:");
                double minIncrease = double.PositiveInfinity;
                for (int i = 0; i < _itemCount; i++)
                {
                    if (!_solution[i] && _weights[i] > 0)
                    {
                        minIncrease = Math.Min(minIncrease, _weights[i] - slack);
                    }
                }
                // Calculate shadow price (marginal value of additional capacity)
                double shadowPrice = 0;
                if (Math.Abs(totalWeight - _capacity) < Epsilon)
                {
                    // If knapsack is full, shadow price is the value/weight ratio of the best excluded item
                    shadowPrice = _values
                        .Select((v, i) => !_solution[i] ? v / _weights[i] : 0)
                        .DefaultIfEmpty(0)
                        .Max();
                    
                    // Store for visualization
                    _shadowPrices["Capacity"] = shadowPrice;
                }  
                if (double.IsFinite(minIncrease) && minIncrease > 0)
                {
                    Console.WriteLine($"  Need to increase capacity by at least {minIncrease:F6} to add another item");
                }
            }
            
            Console.WriteLine($"\nMinimum capacity reduction before solution changes: {minReduction:F6}");
            
            // Calculate criticality ratio for capacity
            if (minReduction > 0 && !double.IsInfinity(minReduction))
            {
                double criticality = minReduction / _capacity * 100;
                Console.WriteLine($"Criticality: The solution remains optimal for a {criticality:F2}% reduction in capacity");
            }

            // Calculate maximum capacity increase without changing the solution
            double maxIncrease = double.PositiveInfinity;
            for (int i = 0; i < _itemCount; i++)
            {
                if (!_solution[i] && _weights[i] > 0)
                {
                    maxIncrease = Math.Min(maxIncrease, _weights[i] - slack);
                }
            }
            if (double.IsPositiveInfinity(maxIncrease))
            {
                Console.WriteLine("No items can be added even with increased capacity");
            }
            else
            {
                Console.WriteLine($"Maximum capacity increase without solution change: {Math.Max(0, maxIncrease):F2}");
            }
        }

        private void AnalyzeValueSensitivity()
        {
            Console.WriteLine("\nVALUE SENSITIVITY");
            Console.WriteLine(new string('-', 40));

            // First, calculate the value-to-weight ratios for all items
            double[] valueToWeight = new double[_itemCount];
            for (int i = 0; i < _itemCount; i++)
            {
                valueToWeight[i] = SafeDivide(_values[i], _weights[i], 0);
            }

            // For each item in the solution, calculate how much its value can change
            for (int i = 0; i < _itemCount; i++)
            {
                if (_solution[i])
                {
                    // For included items, find the maximum possible decrease in value that keeps the item in the solution
                    double maxDecrease = double.PositiveInfinity;
                    
                    // Check against all excluded items that could potentially replace this item
                    for (int j = 0; j < _itemCount; j++)
                    {
                        if (!_solution[j] && _weights[j] <= _capacity - CalculateTotalWeight(_solution, _weights) + _weights[i] + Epsilon)
                        {
                            // Calculate the break-even point where the excluded item becomes more attractive
                            double decrease = _values[i] - (valueToWeight[j] * _weights[i]);
                            if (decrease > 0 && decrease < maxDecrease)
                            {
                                maxDecrease = decrease;
                            }
                        }
                    }
                    
                    string stability = maxDecrease < Epsilon * 10 ? " (unstable)" : "";
                    Console.WriteLine($"Item {i + 1}: Value can decrease by {maxDecrease:F6} before being excluded{stability}");
                }
                else
                {
                    // For excluded items, find the minimum increase needed to include the item
                    double minIncrease = double.PositiveInfinity;
                    
                    // Find the best item to potentially remove to include this one
                    for (int j = 0; j < _itemCount; j++)
                    {
                        if (_solution[j] && _weights[i] <= _capacity - CalculateTotalWeight(_solution, _weights) + _weights[j] + Epsilon)
                        {
                            // Calculate how much we need to increase this item's value to make it worth including
                            double increase = (valueToWeight[j] * _weights[i]) - _values[i];
                            if (increase > 0 && increase < minIncrease)
                            {
                                minIncrease = increase;
                            }
                        }
                    }
                    
                    string potential = double.IsFinite(minIncrease) ? 
                        $" by {minIncrease:F6} to be included" : 
                        " (cannot be included without removing multiple items)";
                    Console.WriteLine($"Item {i + 1}: Value must increase{potential}");
                }
            }
        }

        private void AnalyzeWeightSensitivity()
        {
            Console.WriteLine("\nWEIGHT SENSITIVITY");
            Console.WriteLine(new string('-', 40));

            // Calculate total weight of included items
            double totalWeight = _solution.Select((included, i) => included ? _weights[i] : 0).Sum();
            double remainingCapacity = _capacity - totalWeight;

            for (int i = 0; i < _itemCount; i++)
            {
                string varName = $"Item {i + 1}";
                
                if (_solution[i])
                {
                    // For included items, how much can we decrease the value before it's no longer optimal to include it?
                    double minValue = 0; // Can't go below zero
                    double maxValue = double.PositiveInfinity; // Can increase without bound
                    
                    // Find the best excluded item that could replace this one
                    for (int j = 0; j < _itemCount; j++)
                    {
                        if (!_solution[j] && _weights[j] <= _weights[i])
                        {
                            double ratio = _values[j] / _weights[j];
                            double threshold = ratio * _weights[i];
                            minValue = Math.Max(minValue, threshold);
                        }
                    }
                    
                    // Update objective range for visualization
                    _objectiveRanges[varName] = (minValue, _values[i], maxValue);
                    _reducedCosts[varName] = 0; // Included items have zero reduced cost
                    
                    Console.WriteLine($"Item {i + 1} (included): Value can decrease to {minValue:F2} before being excluded");
                }
                else
                {
                    // For excluded items, how much must the value increase to be included?
                    double threshold = 0;
                    
                    // Find the combination of items that would be replaced
                    double availableCapacity = _capacity - CalculateTotalWeight(_solution, _weights);
                    
                    if (_weights[i] > availableCapacity)
                    {
                        // Need to remove some items to make room
                        // This is a simplification - in practice would need to solve a subproblem
                        threshold = _values[i] / _weights[i];
                    }
                    
                    // Update objective range and reduced cost for visualization
                    _objectiveRanges[varName] = (0, _values[i], threshold);
                    _reducedCosts[varName] = _values[i] - threshold; // Reduced cost is the opportunity cost
                    
                    Console.WriteLine($"Item {i + 1} (excluded): Value must increase to {threshold:F2} to be included");
                }
            }
        }

        private void PrintSubsectionHeader(string title)
        {
            Console.WriteLine("\n" + new string('=', 80));
            Console.WriteLine($"{title}");
            Console.WriteLine(new string('-', title.Length));
        }

        private void AnalyzeCriticalItems()
        {
            PrintSubsectionHeader("CRITICAL ITEMS ANALYSIS");

            var criticalItems = new List<int>();
            var nonCriticalItems = new List<int>();
            
            // Calculate the current total value for comparison
            double currentTotalValue = _solution.Select((included, i) => included ? _values[i] : 0).Sum();
            
            for (int i = 0; i < _itemCount; i++)
            {
                if (!_solution[i]) continue;
                
                // Calculate the best possible solution without this item
                double remainingCapacity = _capacity - CalculateTotalWeight(_solution, _weights) + _weights[i];
                
                // Use a greedy approach to find the best items that fit
                var potentialItems = Enumerable.Range(0, _itemCount)
                    .Where(j => !_solution[j] && j != i) // Exclude current item
                    .Where(j => _weights[j] <= remainingCapacity + Epsilon) // Only items that fit
                    .OrderByDescending(j => _values[j] / _weights[j])
                    .ThenBy(j => _weights[j])
                    .ToList();
                
                double bestValueWithoutI = 0;
                double usedCapacity = 0;
                
                foreach (var j in potentialItems)
                {
                    if (usedCapacity + _weights[j] <= remainingCapacity + Epsilon)
                    {
                        bestValueWithoutI += _values[j];
                        usedCapacity += _weights[j];
                    }
                }
                
                // Calculate the value of the solution without item i
                double solutionWithoutI = currentTotalValue - _values[i] + bestValueWithoutI;
                
                // If the solution without this item is worse, the item is critical
                if (solutionWithoutI < currentTotalValue - Epsilon)
                {
                    criticalItems.Add(i);
                    Console.WriteLine($"Item {i + 1} is CRITICAL: Removing it would reduce value by at least {currentTotalValue - solutionWithoutI:F6}");
                }
                else
                {
                    nonCriticalItems.Add(i);
                    Console.WriteLine($"Item {i + 1} is NOT CRITICAL: Can be removed without loss of value");
                }
            }
            
            // Analyze non-critical items for alternative solutions
            if (nonCriticalItems.Any())
            {
                Console.WriteLine("\nALTERNATIVE SOLUTIONS ANALYSIS");
                Console.WriteLine("The following items can be swapped for others with equal or better value:");
                
                foreach (var i in nonCriticalItems)
                {
                    var alternatives = new List<string>();
                    double remainingCapacity = _capacity - CalculateTotalWeight(_solution, _weights) + _weights[i];
                    
                    // Find items that could replace this one
                    for (int j = 0; j < _itemCount; j++)
                    {
                        if (!_solution[j] && j != i && 
                            _weights[j] <= _weights[i] + remainingCapacity + Epsilon &&
                            _values[j] >= _values[i] - Epsilon)
                        {
                            alternatives.Add($"Item {j + 1} (Value: {_values[j]:F2}, Weight: {_weights[j]:F2})");
                        }
                    }
                    
                    if (alternatives.Any())
                    {
                        Console.WriteLine($"- Item {i + 1} could be replaced with: {string.Join(" or ", alternatives.Take(3))}" +
                                        (alternatives.Count > 3 ? "..." : ""));
                    }
                }
            }
            
            // Final summary
            Console.WriteLine("\nCRITICALITY SUMMARY:");
            Console.WriteLine($"- {criticalItems.Count} critical items (removing any would reduce solution value)");
            Console.WriteLine($"- {nonCriticalItems.Count} non-critical items (can be removed without loss of value)");
            
            if (criticalItems.Count == 0)
            {
                Console.WriteLine("  Note: The solution is highly flexible with multiple optimal configurations.");
            }
            else if (criticalItems.Count == _itemCount)
            {
                Console.WriteLine("  Note: All items are critical - this is a very tight solution.");
            }
        }

        private static double CalculateTotalValue(bool[] solution, double[] values)
        {
            double total = 0;
            for (int i = 0; i < solution.Length; i++)
            {
                if (solution[i]) total += values[i];
            }
            return total;
        }

        private double CalculateTotalWeight(bool[] solution, double[] weights)
        {
            if (solution == null || weights == null || solution.Length != weights.Length)
                return 0;
            double total = 0;
            for (int i = 0; i < solution.Length; i++)
            {
                if (solution[i]) total += weights[i];
            }
            return total;
        }
    }
}
