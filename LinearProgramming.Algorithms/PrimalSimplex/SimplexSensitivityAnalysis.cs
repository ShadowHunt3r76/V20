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
    public class SimplexSensitivityAnalysis : BaseSensitivityAnalysis
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
        public void PerformAnalysis()
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
            int numSlacks = _constraintTypes.Length;
            
            // For each non-basic variable
            for (int j = 0; j < numVars; j++)
            {
                if (!_basisIndices.Contains(j))
                {
                    // Find the column in the tableau corresponding to this variable
                    int col = -1;
                    for (int k = 0; k < _finalTableau.GetLength(1) - 1; k++)
                    {
                        if (_finalTableau[0, k].AlmostEquals(_objectiveCoefficients[j], Epsilon))
                        {
                            col = k;
                            break;
                        }
                    }
                    
                    if (col >= 0)
                    {
                        double reducedCost = _finalTableau[0, col];
                        double currentCoef = _objectiveCoefficients[j];
                        double lowerBound = currentCoef - reducedCost;
                        
                        Console.WriteLine($"{_variableNames[j]}: Current = {currentCoef:F6}");
                        Console.WriteLine($"  Reduced Cost: {reducedCost:F6}");
                        Console.WriteLine($"  Range: [{lowerBound:F6}, +∞)");
                    }
                    else
                    {
                        Console.WriteLine($"{_variableNames[j]}: Could not determine reduced cost - variable not found in tableau");
                    }
                }
            }
            
            // For each basic variable
            for (int i = 0; i < _basisIndices.Length; i++)
            {
                int varIndex = _basisIndices[i];
                if (varIndex < numVars) // Only original variables, not slacks
                {
                    // Calculate range for basic variable's coefficient
                    double lowerBound = double.NegativeInfinity;
                    double upperBound = double.PositiveInfinity;
                    
                    for (int j = 0; j < _finalTableau.GetLength(1) - 1; j++)
                    {
                        if (Math.Abs(_finalTableau[0, j]) > Epsilon && Math.Abs(_finalTableau[i + 1, j]) > Epsilon)
                        {
                            double ratio = _finalTableau[0, j] / _finalTableau[i + 1, j];
                            if (_finalTableau[i + 1, j] > 0)
                            {
                                upperBound = Math.Min(upperBound, ratio);
                            }
                            else
                            {
                                lowerBound = Math.Max(lowerBound, ratio);
                            }
                        }
                    }
                    
                    double currentCoef = _objectiveCoefficients[varIndex];
                    Console.WriteLine($"{_variableNames[varIndex]}: Current = {currentCoef:F6}");
                    Console.WriteLine($"  Range: [{(double.IsNegativeInfinity(lowerBound) ? "-∞" : (currentCoef + lowerBound).ToString("F6"))}, " +
                                    $"{(double.IsPositiveInfinity(upperBound) ? "+∞" : (currentCoef + upperBound).ToString("F6"))}]");
                }
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

        private void AnalyzeReducedCosts()
        {
            Console.WriteLine("\nREDUCED COSTS");
            Console.WriteLine(new string('-', 40));

            int numVars = _variableNames.Length;
            
            for (int j = 0; j < numVars; j++)
            {
                if (!_basisIndices.Contains(j))
                {
                    double reducedCost = 0;
                    // Find the column in the tableau corresponding to this variable
                    for (int k = 0; k < _finalTableau.GetLength(1) - 1; k++)
                    {
                        if (Math.Abs(_finalTableau[0, k] - _objectiveCoefficients[j]) < Epsilon)
                        {
                            reducedCost = _finalTableau[0, k];
                            break;
                        }
                    }
                    
                    Console.WriteLine($"{_variableNames[j]}: {reducedCost:F6}");
                    Console.WriteLine($"  Interpretation: The objective would decrease by {reducedCost:F6} per unit increase in {_variableNames[j]}'s value");
                }
                else
                {
                    Console.WriteLine($"{_variableNames[j]}: 0.000000 (basic variable)");
                }
            }
        }

        private void AnalyzeShadowPrices()
        {
            Console.WriteLine("\nSHADOW PRICES (DUAL VARIABLES)");
            Console.WriteLine(new string('-', 40));

            int numConstraints = _constraintTypes.Length;
            
            for (int i = 0; i < numConstraints; i++)
            {
                double shadowPrice = 0;
                
                // Find the basic variable that corresponds to this constraint's slack
                for (int row = 0; row < _basisIndices.Length; row++)
                {
                    if (_basisIndices[row] == numConstraints + i) // Slack variable for this constraint
                    {
                        shadowPrice = _finalTableau[0, numConstraints + i];
                        break;
                    }
                }
                
                Console.WriteLine($"Constraint {i + 1} ({_constraintTypes[i]}): {shadowPrice:F6}");
                Console.WriteLine($"  Interpretation: The objective would change by {shadowPrice:F6} per unit increase in the RHS");
            }
        }

        private void AnalyzeNonBasicVariables()
        {
            Console.WriteLine("\nNON-BASIC VARIABLES ANALYSIS");
            Console.WriteLine(new string('-', 40));

            int numVars = _variableNames.Length;
            int numConstraints = _constraintTypes.Length;
            
            for (int j = 0; j < numVars; j++)
            {
                if (!_basisIndices.Contains(j))
                {
                    Console.WriteLine($"\n{_variableNames[j]} (non-basic):");
                    
                    // Find the column in the tableau corresponding to this variable
                    for (int k = 0; k < _finalTableau.GetLength(1) - 1; k++)
                    {
                        if (Math.Abs(_finalTableau[0, k] - _objectiveCoefficients[j]) < Epsilon)
                        {
                            double reducedCost = _finalTableau[0, k];
                            Console.WriteLine($"  Reduced Cost: {reducedCost:F6}");
                            
                            // Find the entering variable's impact on basic variables
                            for (int i = 0; i < _basisIndices.Length; i++)
                            {
                                int basisVarIndex = _basisIndices[i];
                                double ratio = _finalTableau[i + 1, k] / _finalTableau[i + 1, _finalTableau.GetLength(1) - 1];
                                
                                if (Math.Abs(_finalTableau[i + 1, k]) > Epsilon)
                                {
                                    string basisVarName = basisVarIndex < _variableNames.Length ? 
                                        _variableNames[basisVarIndex] : $"s{basisVarIndex - numVars + 1}";
                                    
                                    Console.WriteLine($"  Would affect {basisVarName} with ratio {ratio:F6}");
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
}
