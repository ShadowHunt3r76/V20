using System;
using System.Linq;
using Xunit;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.Utils;
using LinearProgramming.Parsing;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Tests.SensitivityAnalysis
{
    public class SimplexSensitivityAnalysisTests : IDisposable
    {
        private readonly double[,] _testTableau;
        private readonly int[] _basisIndices;
        private readonly string[] _variableNames;
        private readonly double[] _objectiveCoefficients;
        private readonly double[] _rhsCoefficients;
        private readonly ConstraintType[] _constraintTypes;

        public SimplexSensitivityAnalysisTests()
        {
            // Standard form of a simple LP problem:
            // Max Z = 3x + 2y
            // Subject to:
            //  2x + y ≤ 100 (slack s1)
            //  x + y ≤ 80   (slack s2)
            //  x ≤ 40       (slack s3)
            //  x, y ≥ 0
            
            // Final optimal tableau (after solving):
            _testTableau = new double[,]
            {
                // x    y    s1   s2   s3   RHS
                { 0.0, 0.0, 1.0, 1.0, 0.0, 20.0 },  // Z equation
                { 0.0, 1.0, 1.0, -1.0, 0.0, 20.0 }, // y equation
                { 0.0, 0.0, -1.0, 2.0, 1.0, 20.0 }, // s3 equation
                { 1.0, 0.0, 0.0, 1.0, 0.0, 40.0 }   // x equation
            };

            // Basis variables in order: Z, y, s3, x
            _basisIndices = new[] { -1, 1, 4, 0 }; // -1 for Z, 1 for y, 4 for s3, 0 for x
            
            _variableNames = new[] { "x", "y", "s1", "s2", "s3" };
            _objectiveCoefficients = new[] { 3.0, 2.0, 0.0, 0.0, 0.0 };
            _rhsCoefficients = new[] { 100.0, 80.0, 40.0 };
            _constraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual };
        }

        [Fact]
        public void Constructor_WithValidParameters_InitializesCorrectly()
        {
            // Act
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
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
            
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            analysis.PerformAnalysis();
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("SIMPLEX SENSITIVITY ANALYSIS", output);
            Assert.Contains("OBJECTIVE COEFFICIENT RANGES", output);
            Assert.Contains("RIGHT-HAND SIDE RANGES", output);
        }

        [Fact]
        public void GetObjectiveCoefficientRange_ForBasicVariable_ReturnsValidRange()
        {
            // Arrange
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act - x is in the basis (index 0)
            var range = analysis.GetObjectiveCoefficientRange(0);

            // Assert - x's coefficient can vary between 2 and 4
            // This is based on the final tableau and reduced costs
            Assert.True(range.lowerBound < 3.0);
            Assert.True(range.upperBound > 3.0);
        }

        [Fact]
        public void GetObjectiveCoefficientRange_ForNonBasicVariable_ReturnsValidRange()
        {
            // Arrange
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act - s1 is not in the basis
            var range = analysis.GetObjectiveCoefficientRange(2);

            // Assert - s1's coefficient must remain <= 1 to maintain optimality
            // (based on the reduced cost of s1 in the final tableau)
            Assert.True(range.lowerBound < 0);
            Assert.True(range.upperBound > 0);
        }

        [Fact]
        public void GetRightHandSideRange_ForActiveConstraint_ReturnsValidRange()
        {
            // Arrange
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
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
        public void IsOptimal_WithValidSolution_ReturnsTrue()
        {
            // Arrange
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
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
        public void GetDualValues_WithValidSolution_ReturnsCorrectDuals()
        {
            // Arrange
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            var duals = analysis.GetDualValues();

            // Assert - Dual values should match the coefficients in the Z row
            // of the final tableau for the slack variables
            Assert.Equal(1.0, duals[0], 6);  // Dual for first constraint
            Assert.Equal(1.0, duals[1], 6);  // Dual for second constraint
            Assert.Equal(0.0, duals[2], 6);  // Dual for third constraint (non-binding)
        }

        [Fact]
        public void GetReducedCosts_WithValidSolution_ReturnsCorrectReducedCosts()
        {
            // Arrange
            var analysis = new SimplexSensitivityAnalysis(
                _testTableau,
                _basisIndices,
                _variableNames,
                _objectiveCoefficients,
                _rhsCoefficients,
                _constraintTypes);

            // Act
            var reducedCosts = analysis.GetReducedCosts();

            // Assert - Non-basic variables should have non-zero reduced costs
            // Basic variables should have reduced cost of 0
            Assert.Equal(0.0, reducedCosts[0], 6);  // x is basic
            Assert.Equal(0.0, reducedCosts[1], 6);  // y is basic
            Assert.True(Math.Abs(reducedCosts[2]) > Epsilon);  // s1 is non-basic
            Assert.True(Math.Abs(reducedCosts[3]) > Epsilon);  // s2 is non-basic
            Assert.Equal(0.0, reducedCosts[4], 6);  // s3 is basic
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
