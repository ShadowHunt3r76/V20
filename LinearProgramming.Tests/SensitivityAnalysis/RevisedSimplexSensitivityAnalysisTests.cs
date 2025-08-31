using System;
using System.IO;
using Xunit;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.Utils;
using LinearProgramming.Parsing;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Tests.SensitivityAnalysis
{
    public class RevisedSimplexSensitivityAnalysisTests : IDisposable
    {
        private readonly double[,] _basisInverse;
        private readonly double[] _basicSolution;
        private readonly double[] _reducedCosts;
        private readonly int[] _basicVariables;
        private readonly string[] _variableNames;
        private readonly double[] _objectiveCoefficients;
        private readonly double[] _rhsCoefficients;
        private readonly ConstraintType[] _constraintTypes;

        public RevisedSimplexSensitivityAnalysisTests()
        {
            // Standard form of a simple LP problem (same as in Simplex test but using revised simplex data structures):
            // Max Z = 3x + 2y
            // Subject to:
            //  2x + y + s1 = 100
            //   x + y + s2 = 80
            //   x + s3 = 40
            //  x, y, s1, s2, s3 ≥ 0
            
            // Optimal basis inverse matrix for the problem
            _basisInverse = new double[,]
            {
                { 1.0, 0.0, 0.0 },
                { 0.0, 1.0, 0.0 },
                { 0.0, 0.0, 1.0 }
            };

            // Basic solution [Z, y, s3, x] (Z is always first in basis)
            _basicSolution = new double[] { 180.0, 20.0, 20.0, 40.0 };
            
            // Reduced costs for [x, y, s1, s2, s3]
            _reducedCosts = new double[] { 0.0, 0.0, 1.0, 1.0, 0.0 };
            
            // Indices of basic variables in order [Z, y, s3, x]
            _basicVariables = new[] { -1, 1, 4, 0 }; // -1 for Z, 1 for y, 4 for s3, 0 for x
            
            _variableNames = new[] { "x", "y", "s1", "s2", "s3" };
            _objectiveCoefficients = new[] { 3.0, 2.0, 0.0, 0.0, 0.0 };
            _rhsCoefficients = new[] { 100.0, 80.0, 40.0 };
            _constraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual };
        }

        [Fact]
        public void Constructor_WithValidParameters_InitializesCorrectly()
        {
            // Act
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                _reducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Assert
            Assert.NotNull(analysis);
        }

        [Fact]
        public void PerformAnalysis_WithValidSolution_OutputsSensitivityInfo()
        {
            // Arrange
            var consoleOutput = new StringWriter();
            Console.SetOut(consoleOutput);
            
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                _reducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            analysis.PerformAnalysis();
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("REVISED SIMPLEX SENSITIVITY ANALYSIS", output);
            Assert.Contains("OBJECTIVE COEFFICIENT RANGES", output);
            Assert.Contains("RIGHT-HAND SIDE RANGES", output);
        }

        [Fact]
        public void GetObjectiveCoefficientRange_ForBasicVariable_ReturnsValidRange()
        {
            // Arrange
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                _reducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act - x is in the basis (index 0)
            var range = analysis.GetObjectiveCoefficientRange(0);

            // Assert - x's coefficient can vary while maintaining optimality
            Assert.True(range.lowerBound < 3.0);
            Assert.True(range.upperBound > 3.0);
        }

        [Fact]
        public void GetRightHandSideRange_ForActiveConstraint_ReturnsValidRange()
        {
            // Arrange
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                _reducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act - First constraint (2x + y ≤ 100) is binding
            var range = analysis.GetRightHandSideRange(0);

            // Assert - The range should be around the current RHS value
            Assert.True(range.lowerBound < 100);
            Assert.True(range.upperBound > 100);
        }

        [Fact]
        public void GetDualValues_WithValidSolution_ReturnsCorrectDuals()
        {
            // Arrange
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                _reducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            var duals = analysis.GetDualValues();

            // Assert - Dual values should match the expected values from the problem
            Assert.Equal(3, duals.Length);
            // First two constraints should have positive dual values (binding)
            Assert.True(duals[0] > 0);
            Assert.True(duals[1] > 0);
            // Third constraint should have zero dual (non-binding)
            Assert.Equal(0.0, duals[2], 6);
        }

        [Fact]
        public void GetReducedCosts_WithValidSolution_ReturnsCorrectReducedCosts()
        {
            // Arrange
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                _reducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            var reducedCosts = analysis.GetReducedCosts();

            // Assert - Should match the input reduced costs
            Assert.Equal(_reducedCosts.Length, reducedCosts.Length);
            for (int i = 0; i < _reducedCosts.Length; i++)
            {
                Assert.Equal(_reducedCosts[i], reducedCosts[i], 6);
            }
        }

        [Fact]
        public void IsOptimal_WithValidSolution_ReturnsTrue()
        {
            // Arrange
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                _reducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            var isOptimal = analysis.IsOptimal();

            // Assert
            Assert.True(isOptimal);
        }

        [Fact]
        public void IsOptimal_WithNegativeReducedCost_ReturnsFalse()
        {
            // Arrange - Create a case where a non-basic variable has negative reduced cost
            var invalidReducedCosts = new[] { 0.0, 0.0, -1.0, 1.0, 0.0 }; // s1 has negative reduced cost
            
            var analysis = new RevisedSimplexSensitivityAnalysis(
                _basisInverse,
                _basicSolution,
                invalidReducedCosts,
                _basicVariables,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            var isOptimal = analysis.IsOptimal();

            // Assert
            Assert.False(isOptimal);
        }

        public void Dispose()
        {
            // Reset console output
            var standardOutput = new StreamWriter(Console.OpenStandardOutput());
            standardOutput.AutoFlush = true;
            Console.SetOut(standardOutput);
        }
    }
}
