using System;

namespace LinearProgramming.Algorithms.Utils
{
    /// <summary>
    /// Additional extension helpers for more readable floating-point comparisons with tolerance.
    /// Consolidates functionality referenced across algorithms, e.g. CompareTo/RoundToEpsilon used in knapsack & cutting-plane code.
    /// </summary>
    public static class MathExtensions
    {
        /// <summary>
        /// Compares two doubles using a tolerance. Behaves like System.Double.CompareTo but treats values within epsilon as equal.
        /// </summary>
        public static int CompareTo(this double value, double other, double epsilon = NumericalStabilityUtils.Epsilon)
        {
            if (Math.Abs(value - other) < epsilon) return 0;
            return value < other ? -1 : 1;
        }

        /// <summary>
        /// Rounds a value to the nearest epsilon grid (effectively zero if within epsilon).
        /// This is useful for avoiding tiny accumulating errors.
        /// </summary>
        public static double RoundToEpsilon(this double value, double epsilon = NumericalStabilityUtils.Epsilon)
        {
            return Math.Abs(value) < epsilon ? 0 : value;
        }
    }
}
