using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.SensitivityAnalysis
{
    /// <summary>
    /// Base class for sensitivity analysis implementations
    /// </summary>
    public abstract class BaseSensitivityAnalysis : ISensitivityAnalysis
    {
        protected const double Epsilon = NumericalStabilityUtils.Epsilon;
        protected readonly string[] _variableNames;
        protected readonly ConstraintType[] _constraintTypes;
        protected readonly int _numVariables;
        protected readonly int _numConstraints;

        protected BaseSensitivityAnalysis(
            string[] variableNames,
            ConstraintType[] constraintTypes,
            int numVariables,
            int numConstraints)
        {
            _variableNames = variableNames ?? throw new ArgumentNullException(nameof(variableNames));
            _constraintTypes = constraintTypes ?? throw new ArgumentNullException(nameof(constraintTypes));
            _numVariables = numVariables;
            _numConstraints = numConstraints;
        }

        protected static string FormatRange(double lower, double upper)
        {
            if (double.IsNegativeInfinity(lower) && double.IsPositiveInfinity(upper))
                return "(-∞, +∞)";
            if (double.IsNegativeInfinity(lower))
                return $"(-∞, {upper:F6}]";
            if (double.IsPositiveInfinity(upper))
                return $"[{lower:F6}, +∞)";
            return $"[{lower:F6}, {upper:F6}]";
        }

        protected static string FormatValueChange(double current, double change, bool isIncrease)
        {
            double newValue = isIncrease ? current + change : current - change;
            string direction = isIncrease ? "increased" : "decreased";
            return $"{current:F6} {direction} by {Math.Abs(change):F6} to {newValue:F6}";
        }

        protected static bool IsWithinTolerance(double value, double target, double tolerance = 0)
        {
            return Math.Abs(value - target) <= (tolerance > 0 ? tolerance : Epsilon);
        }

        protected static bool IsFeasible(double value, double lower, double upper)
        {
            return value >= lower - Epsilon && value <= upper + Epsilon;
        }

        // ISensitivityAnalysis implementation with default implementations
        public abstract void PerformAnalysis();
        
        public virtual (double lowerBound, double upperBound) GetRHSRange(int constraintIndex)
        {
            throw new NotImplementedException("RHS range analysis not implemented for this solver");
        }
        
        public virtual (double lowerBound, double upperBound) GetObjectiveCoefficientRange(int variableIndex)
        {
            throw new NotImplementedException("Objective coefficient range analysis not implemented for this solver");
        }
        
        public virtual double GetReducedCost(int variableIndex)
        {
            throw new NotImplementedException("Reduced cost analysis not implemented for this solver");
        }
        
        public virtual double GetShadowPrice(int constraintIndex)
        {
            throw new NotImplementedException("Shadow price analysis not implemented for this solver");
        }
        
        public virtual Dictionary<string, double> GetSolutionValues()
        {
            throw new NotImplementedException("Solution values not implemented for this solver");
        }
        
        public virtual string GenerateVisualization()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Sensitivity Analysis Results");
            sb.AppendLine(new string('=', 80));
            
            try
            {
                // Basic solution info
                var solution = GetSolutionValues();
                if (solution != null && solution.Count > 0)
                {
                    sb.AppendLine("Solution Values:");
                    foreach (var kvp in solution)
                    {
                        sb.AppendLine($"{kvp.Key}: {kvp.Value:F6}");
                    }
                    sb.AppendLine();
                }
                
                // Objective coefficient ranges
                sb.AppendLine("Objective Coefficient Ranges:");
                for (int i = 0; i < _numVariables && i < _variableNames.Length; i++)
                {
                    try
                    {
                        var range = GetObjectiveCoefficientRange(i);
                        sb.AppendLine($"{_variableNames[i]}: {FormatRange(range.lowerBound, range.upperBound)}");
                    }
                    catch (NotImplementedException)
                    {
                        // Skip if not implemented for this variable
                    }
                }
                
                // RHS ranges
                sb.AppendLine("\nRHS Ranges:");
                for (int i = 0; i < _numConstraints && i < _constraintTypes.Length; i++)
                {
                    try
                    {
                        var range = GetRHSRange(i);
                        sb.AppendLine($"Constraint {i + 1} ({_constraintTypes[i]}): {FormatRange(range.lowerBound, range.upperBound)}");
                    }
                    catch (NotImplementedException)
                    {
                        // Skip if not implemented for this constraint
                    }
                }
                
                // Shadow prices
                sb.AppendLine("\nShadow Prices:");
                for (int i = 0; i < _numConstraints && i < _constraintTypes.Length; i++)
                {
                    try
                    {
                        double price = GetShadowPrice(i);
                        sb.AppendLine($"Constraint {i + 1}: {price:F6}");
                    }
                    catch (NotImplementedException)
                    {
                        // Skip if not implemented
                    }
                }
            }
            catch (Exception ex)
            {
                sb.AppendLine($"Error generating visualization: {ex.Message}");
            }
            
            return sb.ToString();
        }
    }
}
