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
    public class RevisedSimplexSensitivityAnalysis : BaseSensitivityAnalysis
    {
        private readonly double[,] _basisInverse;
        private readonly double[] _reducedCosts;
        private readonly double[] _simplexMultipliers;
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
            _objectiveCoefficients = objectiveCoefficients ?? throw new ArgumentNullException(nameof(objectiveCoefficients));
            _rhsCoefficients = rhsCoefficients ?? throw new ArgumentNullException(nameof(rhsCoefficients));
            _basisIndices = basisIndices ?? throw new ArgumentNullException(nameof(basisIndices));
            _originalVarCount = originalVarCount;

            if (variableNames.Length != objectiveCoefficients.Length)
                throw new ArgumentException("Number of variable names must match number of objective coefficients");
            if (rhsCoefficients.Length != constraintTypes.Length)
                throw new ArgumentException("Number of RHS coefficients must match number of constraints");
        }

        /// <summary>
        /// Performs complete sensitivity analysis on the revised simplex solution
        /// </summary>
        public void PerformAnalysis()
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
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\n⚠ Error during sensitivity analysis: {ex.Message}");
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
            
            for (int i = 0; i < _simplexMultipliers.Length; i++)
            {
                Console.WriteLine($"Constraint {i + 1} ({_constraintTypes[i]}): {_simplexMultipliers[i]:F6}");
                Console.WriteLine($"  Interpretation: The objective would change by {_simplexMultipliers[i]:F6} per unit increase in the RHS");
            }
        }

        private void AnalyzeReducedCosts()
        {
            Console.WriteLine("\nREDUCED COSTS");
            Console.WriteLine(new string('-', 40));

            for (int i = 0; i < _reducedCosts.Length; i++)
            {
                if (i < _variableNames.Length) // Only show original variables
                {
                    string varName = _variableNames[i];
                    double reducedCost = _reducedCosts[i];
                    
                    Console.WriteLine($"{varName}: {reducedCost:F6}");
                    if (Math.Abs(reducedCost) > Epsilon)
                    {
                        Console.WriteLine($"  Interpretation: The objective would change by {reducedCost:F6} per unit increase in {varName}'s value");
                    }
                    else
                    {
                        Console.WriteLine("  This variable is in the optimal basis");
                    }
                }
            }
        }

        private void AnalyzeRHSRanges()
        {
            Console.WriteLine("\nRIGHT-HAND SIDE RANGES");
            Console.WriteLine(new string('-', 40));

            int numConstraints = _constraintTypes.Length;
            
            // Check condition number of basis inverse for numerical stability
            double conditionNumber = NumericalStabilityUtils.ConditionNumber(_basisInverse);
            if (conditionNumber > 1e10)
            {
                Console.WriteLine("⚠ Warning: The basis matrix is ill-conditioned (condition number: {0:E2}). " +
                               "Sensitivity results may be unreliable.", conditionNumber);
            }
            
            for (int i = 0; i < numConstraints; i++)
            {
                double currentRHS = _rhsCoefficients[i];
                double shadowPrice = _simplexMultipliers[i];
                
                // Calculate allowable increase and decrease
                double maxIncrease = double.PositiveInfinity;
                double maxDecrease = double.PositiveInfinity;
                bool hasValidRatios = false;
                
                // For each basic variable, check the ratio with the basis inverse
                for (int j = 0; j < _basisInverse.GetLength(0); j++)
                {
                    double invElement = _basisInverse[j, i];
                    if (!invElement.IsZero(Epsilon))
                    {
                        double ratio = NumericalStabilityUtils.SafeDivide(_currentSolution[j], invElement, Epsilon);
                        if (invElement > Epsilon)
                        {
                            maxDecrease = Math.Min(maxDecrease, ratio);
                            hasValidRatios = true;
                        }
                        else if (invElement < -Epsilon)
                        {
                            maxIncrease = Math.Min(maxIncrease, -ratio);
                            hasValidRatios = true;
                        }
                    }
                }
                
                if (!hasValidRatios)
                {
                    Console.WriteLine($"Constraint {i + 1}: Could not determine valid RHS ranges - numerical instability detected");
                    continue;
                }
                
                Console.WriteLine($"Constraint {i + 1} ({_constraintTypes[i]}):");
                Console.WriteLine($"  Current RHS: {currentRHS:F6}");
                Console.WriteLine($"  Shadow Price: {shadowPrice:F6}");
                
                if (double.IsInfinity(maxDecrease) && double.IsInfinity(maxIncrease))
                {
                    Console.WriteLine("  Range: (-∞, +∞)");
                }
                else
                {
                    double lowerBound = double.IsInfinity(maxDecrease) ? 
                        double.NegativeInfinity : currentRHS - maxDecrease;
                    double upperBound = double.IsInfinity(maxIncrease) ? 
                        double.PositiveInfinity : currentRHS + maxIncrease;
                        
                    Console.WriteLine($"  Range: [{(double.IsNegativeInfinity(lowerBound) ? "-∞" : lowerBound.ToString("F6"))}, " +
                                    $"{(double.IsPositiveInfinity(upperBound) ? "+∞" : upperBound.ToString("F6"))}]");
                }
                
                Console.WriteLine($"  Current Basis Remains Optimal For: " +
                                $"RHS in [{currentRHS - (double.IsInfinity(maxDecrease) ? currentRHS * 0.1 : maxDecrease):F6}, " +
                                $"{currentRHS + (double.IsInfinity(maxIncrease) ? currentRHS * 0.1 : maxIncrease):F6}]");
            }
        }

        private void AnalyzeObjectiveCoefficientRanges()
        {
            Console.WriteLine("\nOBJECTIVE COEFFICIENT RANGES");
            Console.WriteLine(new string('-', 40));

            // For non-basic variables
            for (int j = 0; j < _reducedCosts.Length; j++)
            {
                if (j >= _variableNames.Length) continue; // Skip slack/surplus variables
                
                if (!_basisIndices.Contains(j))
                {
                    double reducedCost = _reducedCosts[j];
                    double currentCoef = _objectiveCoefficients[j];
                    double lowerBound = currentCoef - reducedCost;
                    
                    Console.WriteLine($"{_variableNames[j]} (non-basic):");
                    Console.WriteLine($"  Current Coefficient: {currentCoef:F6}");
                    Console.WriteLine($"  Reduced Cost: {reducedCost:F6}");
                    Console.WriteLine($"  Range: [{lowerBound:F6}, +∞)");
                    Console.WriteLine($"  Interpretation: {_variableNames[j]} will not enter the basis unless its coefficient increases by at least {reducedCost:F6}");
                }
            }
            
            // For basic variables
            Console.WriteLine("\nBasic Variables:");
            for (int i = 0; i < _basisIndices.Length; i++)
            {
                int varIndex = _basisIndices[i];
                if (varIndex < _variableNames.Length) // Only original variables
                {
                    double currentCoef = _objectiveCoefficients[varIndex];
                    double lowerBound = double.NegativeInfinity;
                    double upperBound = double.PositiveInfinity;
                    
                    // Calculate ranges using the basis inverse
                    for (int j = 0; j < _reducedCosts.Length; j++)
                    {
                        if (j >= _variableNames.Length) continue; // Skip slack/surplus variables
                        
                        if (!_basisIndices.Contains(j))
                        {
                            double ratio = _reducedCosts[j] / _basisInverse[i, j];
                            if (_basisInverse[i, j] > 0)
                            {
                                upperBound = Math.Min(upperBound, ratio);
                            }
                            else
                            {
                                lowerBound = Math.Max(lowerBound, ratio);
                            }
                        }
                    }
                    
                    Console.WriteLine($"{_variableNames[varIndex]} (basic):");
                    Console.WriteLine($"  Current Coefficient: {currentCoef:F6}");
                    Console.WriteLine($"  Range: [{(double.IsNegativeInfinity(lowerBound) ? "-∞" : (currentCoef + lowerBound).ToString("F6"))}, " +
                                    $"{(double.IsPositiveInfinity(upperBound) ? "+∞" : (currentCoef + upperBound).ToString("F6"))}]");
                }
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
