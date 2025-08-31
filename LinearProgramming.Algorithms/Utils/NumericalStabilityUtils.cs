using System;
using System.Linq;

namespace LinearProgramming.Algorithms.Utils
{
    /// <summary>
    /// Provides numerical stability utilities for linear programming algorithms.
    /// </summary>
    public static class NumericalStabilityUtils
    {
        /// <summary>
        /// Default epsilon value for floating-point comparisons.
        /// </summary>
        public const double Epsilon = 1e-10;

        /// <summary>
        /// Checks if a value is effectively zero within the specified tolerance.
        /// </summary>
        public static bool IsZero(this double value, double epsilon = Epsilon)
        {
            return Math.Abs(value) < epsilon;
        }

        /// <summary>
        /// Checks if two values are effectively equal within the specified tolerance.
        /// </summary>
        public static bool AlmostEquals(this double a, double b, double epsilon = Epsilon)
        {
            return Math.Abs(a - b) < epsilon;
        }

        /// <summary>
        /// Safely divides two numbers, returning zero if the denominator is effectively zero.
        /// </summary>
        public static double SafeDivide(double numerator, double denominator, double epsilon = Epsilon)
        {
            return denominator.IsZero(epsilon) ? 0 : numerator / denominator;
        }

        /// <summary>
        /// Computes the condition number of a matrix.
        /// </summary>
        public static double ConditionNumber(double[,] matrix)
        {
            // Simple implementation - for more accurate results, consider using a numerical library
            double max = 0;
            double min = double.MaxValue;
            int n = matrix.GetLength(0);
            
            // For simplicity, using the ratio of max to min singular value estimates
            for (int i = 0; i < n; i++)
            {
                double rowSum = 0;
                for (int j = 0; j < n; j++)
                {
                    rowSum += Math.Abs(matrix[i, j]);
                }
                max = Math.Max(max, rowSum);
                min = Math.Min(min, rowSum);
            }
            
            return min > Epsilon ? max / min : double.PositiveInfinity;
        }

        /// <summary>
        /// Checks if a value is within the specified range [lower, upper] with tolerance
        /// </summary>
        public static bool IsInRange(double value, double lower, double upper, double epsilon = Epsilon)
        {
            return value >= lower - epsilon && value <= upper + epsilon;
        }

        /// <summary>
        /// Safely calculates the ratio of two values, handling near-zero cases
        /// </summary>
        public static double SafeRatio(double numerator, double denominator, double defaultValue = 0, double epsilon = Epsilon)
        {
            return Math.Abs(denominator) > epsilon ? numerator / denominator : defaultValue;
        }

        /// <summary>
        /// Determines if two values are approximately equal within the specified tolerance
        /// </summary>
        public static bool AreApproximatelyEqual(double a, double b, double epsilon = Epsilon)
        {
            double absA = Math.Abs(a);
            double absB = Math.Abs(b);
            double diff = Math.Abs(a - b);

            if (a == b)
            {
                return true;
            }
            else if (a == 0 || b == 0 || diff < double.Epsilon)
            {
                return diff < epsilon;
            }
            else
            {
                return diff / (absA + absB) < epsilon;
            }
        }

        /// <summary>
        /// Clamps a value between a minimum and maximum value with optional tolerance
        /// </summary>
        public static double Clamp(double value, double min, double max, double epsilon = Epsilon)
        {
            if (min > max + epsilon)
                throw new ArgumentException($"min ({min}) must be less than or equal to max ({max})");

            if (value < min - epsilon) return min;
            if (value > max + epsilon) return max;
            return value;
        }

        /// <summary>
        /// Rounds a value to zero if it's within epsilon of zero
        /// </summary>
        public static double RoundToZero(double value, double epsilon = Epsilon)
        {
            return Math.Abs(value) < epsilon ? 0 : value;
        }
    }
}
