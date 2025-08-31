using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.SensitivityAnalysis
{
    /// <summary>
    /// Base class for sensitivity analysis implementations
    /// </summary>
    public abstract class BaseSensitivityAnalysis
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
    }
}
