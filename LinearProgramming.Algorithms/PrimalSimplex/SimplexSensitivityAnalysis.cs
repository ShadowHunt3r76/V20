using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Provides sensitivity analysis for Simplex method solutions
    /// </summary>
    public class SimplexSensitivityAnalysis : BaseSensitivityAnalysis, ISensitivityAnalysis
    {
        private readonly double[,] _finalTableau;
        private readonly int[] _basisIndices;
        private readonly double[] _objectiveCoefficients;
        private readonly double[] _rhsCoefficients;

        /// <summary>
        /// Initializes a new instance of the SimplexSensitivityAnalysis class
        /// </summary>
        /// <param name="finalTableau">The final simplex tableau</param>
        /// <param name="basisIndices">Indices of the basis variables</param>
        /// <param name="variableNames">Names of the variables</param>
        /// <param name="objectiveCoefficients">Original objective coefficients</param>
        /// <param name="rhsCoefficients">Right-hand side coefficients</param>
        /// <param name="constraintTypes">Constraint types (≤, =, ≥)</param>
        public SimplexSensitivityAnalysis(
            double[,] finalTableau,
            int[] basisIndices,
            string[] variableNames,
            double[] objectiveCoefficients,
            double[] rhsCoefficients,
            ConstraintType[] constraintTypes)
            : base(variableNames, constraintTypes, variableNames?.Length ?? 0, constraintTypes?.Length ?? 0)
        {
            _finalTableau = finalTableau ?? throw new ArgumentNullException(nameof(finalTableau));
            _basisIndices = basisIndices ?? throw new ArgumentNullException(nameof(basisIndices));
            _objectiveCoefficients = objectiveCoefficients ?? throw new ArgumentNullException(nameof(objectiveCoefficients));
            _rhsCoefficients = rhsCoefficients ?? throw new ArgumentNullException(nameof(rhsCoefficients));

            if (variableNames.Length != objectiveCoefficients.Length)
                throw new ArgumentException("Number of variable names must match number of objective coefficients");
            if (rhsCoefficients.Length != constraintTypes.Length)
                throw new ArgumentException("Number of RHS coefficients must match number of constraints");
        }

        /// <summary>
        /// Performs complete sensitivity analysis on the simplex solution
        /// </summary>
        public override void PerformAnalysis()
        {
            Console.WriteLine("\n" + new string('=', 80));
            Console.WriteLine("SIMPLEX SENSITIVITY ANALYSIS");
            Console.WriteLine(new string('=', 80));

            try
            {
                // 1. Basic solution information
                DisplaySolutionInfo();

                // 2. Objective function coefficient ranges
                AnalyzeObjectiveCoefficientRanges();

                // 3. Right-hand side ranges
                AnalyzeRHSRanges();

                // 4. Reduced costs analysis
                AnalyzeReducedCosts();

                // 5. Shadow prices (dual variables)
                AnalyzeShadowPrices();

                // 6. Allowable ranges for non-basic variables
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
            Console.WriteLine("\nBasic Variables:");
            for (int i = 0; i < _basisIndices.Length; i++)
            {
                int varIndex = _basisIndices[i];
                string varName = varIndex < _variableNames.Length ? _variableNames[varIndex] : $"x{varIndex + 1}";
                double value = _finalTableau[i + 1, _finalTableau.GetLength(1) - 1];
                Console.WriteLine($"{varName} = {value:F6}");
            }
            
            // Objective value
            double objValue = _finalTableau[0, _finalTableau.GetLength(1) - 1];
            Console.WriteLine($"\nOptimal Objective Value: {objValue:F6}");
        }

        private void AnalyzeObjectiveCoefficientRanges()
        {
            Console.WriteLine("\nOBJECTIVE COEFFICIENT RANGES");
            Console.WriteLine(new string('-', 40));

            int numVars = _variableNames.Length;
            var objectiveRanges = new Dictionary<string, (double min, double current, double max)>();
            
            // For each variable
            for (int varIndex = 0; varIndex < numVars; varIndex++)
            {
                string varName = _variableNames[varIndex];
                double currentCoefficient = _objectiveCoefficients[varIndex];
                
                if (_basisIndices.Contains(varIndex))
                {
                    // Basic variable - calculate allowable increase and decrease
                    double maxIncrease = double.PositiveInfinity;
                    double maxDecrease = double.PositiveInfinity;
                    int basisRow = Array.IndexOf(_basisIndices, varIndex) + 1; // +1 because of objective row
                    
                    // For each non-basic variable
                    for (int j = 0; j < _finalTableau.GetLength(1) - 1; j++)
                    {
                        if (!_basisIndices.Contains(j) && !_finalTableau[basisRow, j].AlmostEquals(0, Epsilon))
                        {
                            double ratio = -_finalTableau[0, j] / _finalTableau[basisRow, j];
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
                    objectiveRanges[varName] = (
                        min: currentCoefficient - maxDecrease,
                        current: currentCoefficient,
                        max: currentCoefficient + maxIncrease
                    );
                }
                else
                {
                    // Non-basic variable - range is from -∞ to current + reduced cost
                    double reducedCost = _finalTableau[0, varIndex];
                    double upperBound = currentCoefficient + reducedCost;
                    
                    Console.WriteLine($"{varName} (non-basic):");
                    Console.WriteLine($"  Current coefficient: {currentCoefficient:F6}");
                    Console.WriteLine($"  Reduced cost: {reducedCost:F6}");
                    Console.WriteLine($"  Range: (-∞, {upperBound:F6}]");
                    
                    // Store for visualization
                    objectiveRanges[varName] = (
                        min: double.NegativeInfinity,
                        current: currentCoefficient,
                        max: upperBound
                    );
                }
            }
            
            // Add visualization for objective coefficient ranges
            if (objectiveRanges.Count > 0)
            {
                Console.WriteLine("\n" + SensitivityVisualizer.CreateBarChart(
                    "OBJECTIVE COEFFICIENT RANGES",
                    objectiveRanges,
                    width: 50
                ));
            }
        }

        private void AnalyzeRHSRanges()
        {
            Console.WriteLine("\nRIGHT-HAND SIDE RANGES");
            Console.WriteLine(new string('-', 40));

            int numConstraints = _constraintTypes.Length;
            int rhsCol = _finalTableau.GetLength(1) - 1;
            
            for (int i = 0; i < numConstraints; i++)
            {
                double currentRHS = _rhsCoefficients[i];
                double shadowPrice = 0;
                double lowerBound = double.NegativeInfinity;
                double upperBound = double.PositiveInfinity;
                
                // Find the basic variable in this row
                int basisRow = -1;
                for (int row = 0; row < _basisIndices.Length; row++)
                {
                    if (_basisIndices[row] == numConstraints + i) // Slack variable for this constraint
                    {
                        basisRow = row + 1; // +1 because of objective row
                        shadowPrice = _finalTableau[0, numConstraints + i];
                        break;
                    }
                }
                
                if (basisRow > 0)
                {
                    // Calculate range for RHS
                    for (int col = 0; col < _finalTableau.GetLength(1) - 1; col++)
                    {
                        if (Math.Abs(_finalTableau[basisRow, col]) > Epsilon)
                        {
                            double ratio = _finalTableau[basisRow, rhsCol] / _finalTableau[basisRow, col];
                            if (_finalTableau[basisRow, col] > 0)
                            {
                                upperBound = Math.Min(upperBound, ratio);
                            }
                            else
                            {
                                lowerBound = Math.Max(lowerBound, ratio);
                            }
                        }
                    }
                }
                
                Console.WriteLine($"Constraint {i + 1} ({_constraintTypes[i]}):");
                Console.WriteLine($"  Current RHS: {currentRHS:F6}");
                Console.WriteLine($"  Shadow Price: {shadowPrice:F6}");
                
                if (double.IsNegativeInfinity(lowerBound) && double.IsPositiveInfinity(upperBound))
                {
                    Console.WriteLine("  Range: (-∞, +∞)");
                }
                else
                {
                    Console.WriteLine($"  Range: [{(double.IsNegativeInfinity(lowerBound) ? "-∞" : (currentRHS + lowerBound).ToString("F6"))}, " +
                                    $"{(double.IsPositiveInfinity(upperBound) ? "+∞" : (currentRHS + upperBound).ToString("F6"))}]");
                }
            }
        }

        private Dictionary<string, double> _reducedCosts = new();
        
        private void AnalyzeReducedCosts()
        {
            Console.WriteLine("\nREDUCED COSTS");
            Console.WriteLine(new string('-', 40));
            
            // Clear previous reduced costs
            _reducedCosts.Clear();
            var nonZeroReducedCosts = new Dictionary<string, double>();

            for (int j = 0; j < _variableNames.Length; j++)
            {
                if (!_basisIndices.Contains(j)) // Only for non-basic variables
                {
                    double reducedCost = _finalTableau[0, j];
                    _reducedCosts[_variableNames[j]] = reducedCost;
                    
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

        private Dictionary<string, double> _shadowPrices = new();
        
        private void AnalyzeShadowPrices()
        {
            Console.WriteLine("\nSHADOW PRICES");
            Console.WriteLine(new string('-', 40));

            int numSlacks = _finalTableau.GetLength(1) - _variableNames.Length - 1; // Exclude RHS column and original variables
            
            // Clear previous shadow prices
            _shadowPrices.Clear();
            
            for (int i = 0; i < numSlacks && i < _constraintTypes.Length; i++)
            {
                int slackIndex = _variableNames.Length + i;
                double shadowPrice = _finalTableau[0, slackIndex];
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

        private void AnalyzeNonBasicVariables()
        {
            Console.WriteLine("\nNON-BASIC VARIABLES ANALYSIS");
            Console.WriteLine(new string('-', 40));

            int numVars = _finalTableau.GetLength(1) - 1; // Exclude RHS column
            int numConstraints = _finalTableau.GetLength(0) - 1; // Exclude objective row

            // Prepare data for visualization
            var nonBasicVars = new Dictionary<string, (double min, double current, double max)>();
            
            for (int j = 0; j < numVars; j++)
            {
                if (!_basisIndices.Contains(j))
                {
                    string varName = j < _variableNames.Length ? _variableNames[j] : $"x{j + 1}";
                    double reducedCost = _finalTableau[0, j];
                    
                    Console.WriteLine($"{varName} (reduced cost = {reducedCost:F6}):");
                    Console.WriteLine($"  Current value: 0 (non-basic)");
                    
                    // Calculate allowable increase/decrease for this non-basic variable
                    double maxIncrease = double.PositiveInfinity;
                    double maxDecrease = double.PositiveInfinity;
                    
                    for (int i = 0; i < numConstraints; i++)
                    {
                        double a_ij = _finalTableau[i + 1, j];
                        double b_i = _finalTableau[i + 1, numVars];
                        
                        if (a_ij > Epsilon)
                        {
                            double ratio = b_i / a_ij;
                            maxIncrease = Math.Min(maxIncrease, ratio);
                        }
                        else if (a_ij < -Epsilon)
                        {
                            double ratio = b_i / a_ij;
                            maxDecrease = Math.Min(maxDecrease, -ratio);
                        }
                    }
                    
                    Console.WriteLine($"  Allowable increase: {(double.IsInfinity(maxIncrease) ? "∞" : maxIncrease.ToString("F6"))}");
                    Console.WriteLine($"  Allowable decrease: {(double.IsInfinity(maxDecrease) ? "∞" : maxDecrease.ToString("F6"))}");
                    
                    // Store for visualization
                    nonBasicVars[varName] = (
                        min: -maxDecrease,
                        current: 0,
                        max: maxIncrease
                    );
                }
            }
            
            // Add visualization for non-basic variables
            if (nonBasicVars.Count > 0)
            {
                Console.WriteLine("\n" + SensitivityVisualizer.CreateBarChart(
                    "NON-BASIC VARIABLES RANGES",
                    nonBasicVars,
                    width: 40
                ));
            }
        }
    }
}
