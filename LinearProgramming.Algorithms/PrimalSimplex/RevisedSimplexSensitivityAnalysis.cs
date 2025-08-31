using System;
using System.Linq;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Provides sensitivity analysis for Revised Primal Simplex method solutions
    /// </summary>
    public class RevisedSimplexSensitivityAnalysis : BaseSensitivityAnalysis, ISensitivityAnalysis
    {
        private readonly double[,] _basisInverse;
        private readonly double[] _reducedCosts;
        private readonly double[] _simplexMultipliers;
        private readonly double[] _currentSolution;
        private readonly double[] _objectiveCoefficients;
        private readonly double[] _rhsCoefficients;
        private readonly int[] _basisIndices;
        private readonly int _originalVarCount;

        /// <summary>
        /// Initializes a new instance of the RevisedSimplexSensitivityAnalysis class
        /// </summary>
        /// <param name="basisInverse">Inverse of the final basis matrix</param>
        /// <param name="reducedCosts">Reduced costs from the final solution</param>
        /// <param name="simplexMultipliers">Dual variables (simplex multipliers)</param>
        /// <param name="objectiveCoefficients">Original objective coefficients</param>
        /// <param name="rhsCoefficients">Right-hand side coefficients</param>
        /// <param name="variableNames">Names of the variables</param>
        /// <param name="constraintTypes">Constraint types (≤, =, ≥)</param>
        /// <param name="basisIndices">Indices of the basis variables</param>
        /// <param name="originalVarCount">Number of original variables (before adding slacks/artificials)</param>
        public RevisedSimplexSensitivityAnalysis(
            double[,] basisInverse,
            double[] reducedCosts,
            double[] simplexMultipliers,
            double[] currentSolution,
            double[] objectiveCoefficients,
            double[] rhsCoefficients,
            string[] variableNames,
            ConstraintType[] constraintTypes,
            int[] basisIndices,
            int originalVarCount)
            : base(variableNames, constraintTypes, variableNames?.Length ?? 0, constraintTypes?.Length ?? 0)
        {
            _basisInverse = basisInverse ?? throw new ArgumentNullException(nameof(basisInverse));
            _reducedCosts = reducedCosts ?? throw new ArgumentNullException(nameof(reducedCosts));
            _simplexMultipliers = simplexMultipliers ?? throw new ArgumentNullException(nameof(simplexMultipliers));
            _currentSolution = currentSolution ?? throw new ArgumentNullException(nameof(currentSolution));
            _objectiveCoefficients = objectiveCoefficients ?? throw new ArgumentNullException(nameof(objectiveCoefficients));
            _rhsCoefficients = rhsCoefficients ?? throw new ArgumentNullException(nameof(rhsCoefficients));
            _basisIndices = basisIndices ?? throw new ArgumentNullException(nameof(basisIndices));
            _originalVarCount = originalVarCount;

            if (variableNames.Length != objectiveCoefficients.Length)
                throw new ArgumentException("Number of variable names must match number of objective coefficients");
            if (rhsCoefficients.Length != constraintTypes.Length)
                throw new ArgumentException("Number of RHS coefficients must match number of constraints");
        }

        private readonly Dictionary<string, (double min, double current, double max)> _objectiveRanges = new();
        private readonly Dictionary<string, (double min, double current, double max)> _rhsRanges = new();
        private readonly Dictionary<string, double> _shadowPrices = new();
        private readonly Dictionary<string, double> _reducedCostsDict = new();
        
        /// <summary>
        /// Performs complete sensitivity analysis on the revised simplex solution
        /// </summary>
        public override void PerformAnalysis()
        {
            Console.WriteLine("\n" + new string('=', 80));
            Console.WriteLine("REVISED SIMPLEX SENSITIVITY ANALYSIS");
            Console.WriteLine(new string('=', 80));

            try
            {
                // 1. Basic solution information
                DisplaySolutionInfo();

                // 2. Shadow prices (dual variables)
                AnalyzeShadowPrices();

                // 3. Reduced costs analysis
                AnalyzeReducedCosts();

                // 4. Right-hand side ranges
                AnalyzeRHSRanges();

                // 5. Objective coefficient ranges
                AnalyzeObjectiveCoefficientRanges();

                // 6. Non-basic variables analysis
                AnalyzeNonBasicVariables();
                
                // 7. Generate and display visualization
                Console.WriteLine("\n" + new string('-', 40));
                Console.WriteLine("VISUALIZATION");
                Console.WriteLine(new string('-', 40));
                Console.WriteLine(GenerateVisualization());
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\n⚠ Error during sensitivity analysis: {ex.Message}");
                if (ex.InnerException != null)
                {
                    Console.WriteLine($"Inner exception: {ex.InnerException.Message}");
                }
            }

            Console.WriteLine("\n" + new string('=', 80));
            Console.WriteLine("END OF SENSITIVITY ANALYSIS");
            Console.WriteLine(new string('=', 80));
        }

        private void DisplaySolutionInfo()
        {
            Console.WriteLine("\nSOLUTION INFORMATION");
            Console.WriteLine(new string('-', 40));
            
            // Basic variables and their values
            Console.WriteLine("\nBasic Variables and Values:");
            for (int i = 0; i < _basisIndices.Length; i++)
            {
                int varIndex = _basisIndices[i];
                string varName = varIndex < _variableNames.Length ? _variableNames[varIndex] : $"x{varIndex + 1}";
                Console.WriteLine($"{varName} = {_currentSolution[i]:F6}");
            }
            
            // Objective value
            double objValue = CalculateCurrentObjective();
            Console.WriteLine($"\nOptimal Objective Value: {objValue:F6}");
        }

        private double CalculateCurrentObjective()
        {
            double objValue = 0;
            for (int i = 0; i < _basisIndices.Length; i++)
            {
                int varIndex = _basisIndices[i];
                if (varIndex < _objectiveCoefficients.Length)
                {
                    objValue += _objectiveCoefficients[varIndex] * _currentSolution[i];
                }
            }
            return objValue;
        }

        private void AnalyzeShadowPrices()
        {
            Console.WriteLine("\nSHADOW PRICES (DUAL VARIABLES)");
            Console.WriteLine(new string('-', 40));
            
            // Clear previous shadow prices
            _shadowPrices.Clear();

            for (int i = 0; i < _simplexMultipliers.Length; i++)
            {
                double shadowPrice = _simplexMultipliers[i];
                string constraintName = $"Constraint {i + 1} ({_constraintTypes[i]})";
                
                Console.WriteLine($"{constraintName}: {shadowPrice:F6}");
                Console.WriteLine($"  Interpretation: The objective would change by {shadowPrice:F6} per unit increase in the RHS");
                
                // Store for visualization
                _shadowPrices[constraintName] = shadowPrice;
            }
            
            // Add visualization for shadow prices
            if (_shadowPrices.Count > 0)
            {
                Console.WriteLine("\n" + SensitivityVisualizer.CreateHistogram(
                    "SHADOW PRICES (DUAL VARIABLES)",
                    _shadowPrices,
                    height: 10,
                    width: 50
                ));
            }
        }

        private void AnalyzeReducedCosts()
        {
            Console.WriteLine("\nREDUCED COSTS");
            Console.WriteLine(new string('-', 40));
            
            // Clear previous reduced costs
            _reducedCostsDict.Clear();
            var nonZeroReducedCosts = new Dictionary<string, double>();

            for (int j = 0; j < _originalVarCount; j++)
            {
                if (!_basisIndices.Contains(j)) // Only for non-basic variables
                {
                    double reducedCost = _reducedCosts[j];
                    _reducedCostsDict[_variableNames[j]] = reducedCost;
                    
                    if (Math.Abs(reducedCost) > Epsilon)
                    {
                        Console.WriteLine($"{_variableNames[j]}: {reducedCost:F6}");
                        Console.WriteLine($"  Interpretation: The objective would change by {reducedCost:F6} if this variable enters the basis");
                        nonZeroReducedCosts[_variableNames[j]] = reducedCost;
                    }
                }
            }
            
            // Add visualization for non-zero reduced costs
            if (nonZeroReducedCosts.Count > 0)
            {
                Console.WriteLine("\n" + SensitivityVisualizer.CreateHistogram(
                    "NON-ZERO REDUCED COSTS",
                    nonZeroReducedCosts,
                    height: 10,
                    width: 50
                ));
            }
        }

        private void AnalyzeRHSRanges()
        {
            Console.WriteLine("\nRIGHT-HAND SIDE RANGES");
            Console.WriteLine(new string('-', 40));

            int numConstraints = _constraintTypes.Length;
            
            // Clear previous RHS ranges
            _rhsRanges.Clear();
            
            for (int i = 0; i < numConstraints; i++)
            {
                double currentRHS = _rhsCoefficients[i];
                double shadowPrice = _simplexMultipliers[i];
                string constraintName = $"Constraint {i + 1} ({_constraintTypes[i]})";
                
                Console.WriteLine($"{constraintName}:");
                Console.WriteLine($"  Current RHS: {currentRHS:F6}");
                Console.WriteLine($"  Shadow price: {shadowPrice:F6}");
                
                // Calculate allowable increase/decrease for this RHS
                double maxIncrease = double.PositiveInfinity;
                double maxDecrease = double.PositiveInfinity;
                
                for (int j = 0; j < numConstraints; j++)
                {
                    double a_ij = _basisInverse[j, i];
                    double b_i = _currentSolution[j];
                    
                    if (a_ij > Epsilon)
                    {
                        double ratio = b_i / a_ij;
                        maxDecrease = Math.Min(maxDecrease, ratio);
                    }
                    else if (a_ij < -Epsilon)
                    {
                        double ratio = b_i / a_ij;
                        maxIncrease = Math.Min(maxIncrease, -ratio);
                    }
                }
                
                Console.WriteLine($"  Allowable increase: {(double.IsInfinity(maxIncrease) ? "∞" : maxIncrease.ToString("F6"))}");
                Console.WriteLine($"  Allowable decrease: {(double.IsInfinity(maxDecrease) ? "∞" : maxDecrease.ToString("F6"))}");
                Console.WriteLine($"  Range: [{currentRHS - maxDecrease:F6}, {currentRHS + maxIncrease:F6}]");
                
                // Store for visualization
                _rhsRanges[constraintName] = (
                    min: currentRHS - maxDecrease,
                    current: currentRHS,
                    max: currentRHS + maxIncrease
                );
            }
            
            // Add visualization for RHS ranges
            if (_rhsRanges.Count > 0)
            {
                Console.WriteLine("\n" + SensitivityVisualizer.CreateBarChart(
                    "RIGHT-HAND SIDE RANGES",
                    _rhsRanges,
                    width: 50
                ));
            }
        }

        private void AnalyzeObjectiveCoefficientRanges()
        {
            Console.WriteLine("\nOBJECTIVE COEFFICIENT RANGES");
            Console.WriteLine(new string('-', 40));
            
            // Clear previous objective ranges
            _objectiveRanges.Clear();

            for (int j = 0; j < _originalVarCount; j++)
            {
                string varName = _variableNames[j];
                double currentCoefficient = _objectiveCoefficients[j];
                
                if (_basisIndices.Contains(j))
                {
                    // Basic variable - calculate allowable increase/decrease
                    double maxIncrease = double.PositiveInfinity;
                    double maxDecrease = double.PositiveInfinity;
                    
                    for (int k = 0; k < _originalVarCount; k++)
                    {
                        if (k != j && !_basisIndices.Contains(k))
                        {
                            double ratio = -_reducedCosts[k] / _basisInverse[Array.IndexOf(_basisIndices, j), k];
                            if (ratio > 0)
                            {
                                maxIncrease = Math.Min(maxIncrease, ratio);
                            }
                            else
                            {
                                maxDecrease = Math.Min(maxDecrease, -ratio);
                            }
                        }
                    }
                    
                    Console.WriteLine($"{varName} (basic):");
                    Console.WriteLine($"  Current coefficient: {currentCoefficient:F6}");
                    Console.WriteLine($"  Allowable increase: {(double.IsInfinity(maxIncrease) ? "∞" : maxIncrease.ToString("F6"))}");
                    Console.WriteLine($"  Allowable decrease: {(double.IsInfinity(maxDecrease) ? "∞" : maxDecrease.ToString("F6"))}");
                    Console.WriteLine($"  Range: [{currentCoefficient - maxDecrease:F6}, {currentCoefficient + maxIncrease:F6}]");
                    
                    // Store for visualization
                    _objectiveRanges[varName] = (
                        min: currentCoefficient - maxDecrease,
                        current: currentCoefficient,
                        max: currentCoefficient + maxIncrease
                    );
                }
                else
                {
                    // Non-basic variable - range is from -∞ to current + reduced cost
                    double upperBound = currentCoefficient + _reducedCosts[j];
                    Console.WriteLine($"{varName} (non-basic):");
                    Console.WriteLine($"  Current coefficient: {currentCoefficient:F6}");
                    Console.WriteLine($"  Reduced cost: {_reducedCosts[j]:F6}");
                    Console.WriteLine($"  Range: (-∞, {upperBound:F6}]");
                    
                    // Store for visualization
                    _objectiveRanges[varName] = (
                        min: double.NegativeInfinity,
                        current: currentCoefficient,
                        max: upperBound
                    );
                }
            }
            
            // Add visualization for objective coefficient ranges
            if (_objectiveRanges.Count > 0)
            {
                Console.WriteLine("\n" + SensitivityVisualizer.CreateBarChart(
                    "OBJECTIVE COEFFICIENT RANGES",
                    _objectiveRanges,
                    width: 50
                ));
            }
        }

        private void AnalyzeNonBasicVariables()
        {
            Console.WriteLine("\nNON-BASIC VARIABLES ANALYSIS");
            Console.WriteLine(new string('-', 40));

            for (int j = 0; j < _variableNames.Length; j++)
            {
                if (!_basisIndices.Contains(j))
                {
                    Console.WriteLine($"\n{_variableNames[j]} (non-basic):");
                    double reducedCost = _reducedCosts[j];
                    Console.WriteLine($"  Reduced Cost: {reducedCost:F6}");
                    
                    // Calculate impact on basic variables if this variable enters the basis
                    for (int i = 0; i < _basisIndices.Length; i++)
                    {
                        int basisVarIndex = _basisIndices[i];
                        if (Math.Abs(_basisInverse[i, j]) > Epsilon)
                        {
                            double ratio = _currentSolution[i] / _basisInverse[i, j];
                            string basisVarName = basisVarIndex < _variableNames.Length ? 
                                _variableNames[basisVarIndex] : $"s{basisVarIndex - _originalVarCount + 1}";
                            
                            Console.WriteLine($"  Would affect {basisVarName} with ratio {ratio:F6}");
                        }
                    }
                }
            }
        }
    }
}
