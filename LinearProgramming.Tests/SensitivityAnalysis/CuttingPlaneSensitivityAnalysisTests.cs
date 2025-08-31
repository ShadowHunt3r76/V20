using System;
using System.IO;
using Xunit;
using LinearProgramming.Algorithms.CuttingPlane;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Parsing;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Tests.SensitivityAnalysis
{
    public class CuttingPlaneSensitivityAnalysisTests : IDisposable
    {
        private readonly CuttingPlaneSolution _testSolution;
        private readonly ParsedLinearProgrammingModel _testModel;

        public CuttingPlaneSensitivityAnalysisTests()
        {
            // Create a test model (Integer Programming problem)
            _testModel = new ParsedLinearProgrammingModel
            {
                Objective = new Objective { Type = OptimizationType.Maximize },
                Variables = new[]
                {
                    new Variable { Name = "x1", Type = VariableType.Integer },
                    new Variable { Name = "x2", Type = VariableType.Integer }
                },
                ObjectiveCoefficients = new[] { 3.0, 2.0 },
                ConstraintCoefficients = new[] { new[] { 2.0, 1.0 }, new[] { 1.0, 1.0 } },
                ConstraintConstants = new[] { 10.0, 6.0 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual }
            };

            // Create a test solution
            _testSolution = new CuttingPlaneSolution
            {
                BestObjectiveValue = 12.0,  // x1=2, x2=3
                BestSolution = new[] { 2.0, 3.0 },
                OverallStatus = SolutionStatus.Optimal,
                Iterations = 5,
                CutsAdded = 2,
                FinalLPSolution = new PrimalSimplexSolution
                {
                    Status = SolutionStatus.Optimal,
                    ObjectiveValue = 12.0,
                    VariableValues = new[] { 2.0, 3.0 },
                    // Mock optimal tableau for the final LP relaxation
                    OptimalTable = new double[,]
                    {
                        { 1.0, 0.0, 1.0, 1.0, 12.0 },  // Z + s1 + s2 = 12
                        { 0.0, 1.0, 1.0, -1.0, 2.0 },   // x1 + s1 - s2 = 2
                        { 0.0, 0.0, -1.0, 2.0, 3.0 }    // -s1 + 2s2 = 3 (x2)
                    }
                },
                Model = new ParsedLinearProgrammingModel
                {
                    Variables = _testModel.Variables,
                    ObjectiveCoefficients = _testModel.ObjectiveCoefficients,
                    ConstraintTypes = _testModel.ConstraintTypes,
                    ConstraintCoefficients = _testModel.ConstraintCoefficients,
                    ConstraintConstants = _testModel.ConstraintConstants
                },
                CutHistory = new[]
                {
                    new CutInfo { Iteration = 1, CutCoefficients = new[] { -0.5, 1.0 }, Rhs = 1.5 },
                    new CutInfo { Iteration = 3, CutCoefficients = new[] { 1.0, 0.0 }, Rhs = 2.0 }
                }
            };
        }

        [Fact]
        public void Constructor_WithValidParameters_InitializesCorrectly()
        {
            // Act
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);

            // Assert
            Assert.NotNull(analysis);
        }

        [Fact]
        public void PerformAnalysis_WithValidSolution_OutputsSensitivityInfo()
        {
            // Arrange
            var consoleOutput = new StringWriter();
            Console.SetOut(consoleOutput);
            
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);

            // Act
            analysis.PerformAnalysis();
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("CUTTING PLANE SENSITIVITY ANALYSIS", output);
            Assert.Contains("SOLUTION INFORMATION", output);
            Assert.Contains("LP-BASED SENSITIVITY ANALYSIS", output);
            Assert.Contains("CUT-SPECIFIC SENSITIVITY ANALYSIS", output);
        }

        [Fact]
        public void PerformCutSpecificAnalysis_WithValidSolution_OutputsCutInfo()
        {
            // Arrange
            var consoleOutput = new StringWriter();
            Console.SetOut(consoleOutput);
            
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);
            
            // Use reflection to access the protected method for testing
            var methodInfo = typeof(CuttingPlaneSensitivityAnalysis)
                .GetMethod("PerformCutSpecificAnalysis", 
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            // Act
            methodInfo.Invoke(analysis, null);
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("CUT-SPECIFIC SENSITIVITY ANALYSIS", output);
            Assert.Contains("Iteration", output);
            Assert.Contains("Impact", output);
        }

        [Fact]
        public void GetObjectiveCoefficientRange_WithValidVariable_ReturnsValidRange()
        {
            // Arrange
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);

            // Act - Test with x1 (in the solution)
            var range = analysis.GetObjectiveCoefficientRange(0);

            // Assert - The range should be reasonable values around the original coefficient (3.0)
            Assert.True(range.lowerBound < 3.0);
            Assert.True(range.upperBound > 3.0);
        }

        [Fact]
        public void GetRightHandSideRange_WithValidConstraint_ReturnsValidRange()
        {
            // Arrange
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);

            // Act - Test with the first constraint (2x1 + x2 ≤ 10)
            var range = analysis.GetRightHandSideRange(0);

            // Assert - The range should be reasonable values around the original RHS (10.0)
            Assert.True(range.lowerBound < 10.0);
            Assert.True(range.upperBound > 10.0);
        }

        [Fact]
        public void GetDualValues_WithValidSolution_ReturnsCorrectDuals()
        {
            // Arrange
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var duals = analysis.GetDualValues();

            // Assert - Should return dual values for all original constraints and cuts
            // Original constraints + cuts (2 original + 2 cuts in this case)
            Assert.Equal(4, duals.Length);
            // Duals should be non-negative for ≤ constraints
            Assert.All(duals, d => Assert.True(d >= 0));
        }

        [Fact]
        public void GetReducedCosts_WithValidSolution_ReturnsCorrectReducedCosts()
        {
            // Arrange
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var reducedCosts = analysis.GetReducedCosts();

            // Assert - Should return reduced costs for all variables
            Assert.Equal(2, reducedCosts.Length); // Two variables in our test case
            // Basic variables should have reduced cost of 0
            // In our test case, both x1 and x2 are in the basis
            Assert.Equal(0.0, reducedCosts[0], 6); // x1 is basic
            Assert.Equal(0.0, reducedCosts[1], 6); // x2 is basic
        }

        [Fact]
        public void IsOptimal_WithValidSolution_ReturnsTrue()
        {
            // Arrange
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var isOptimal = analysis.IsOptimal();

            // Assert
            Assert.True(isOptimal);
        }

        [Fact]
        public void GetCutImpact_WithValidCut_ReturnsCorrectImpact()
        {
            // Arrange
            var analysis = new CuttingPlaneSensitivityAnalysis(_testSolution, _testModel);
            
            // Use reflection to access the protected method for testing
            var methodInfo = typeof(CuttingPlaneSensitivityAnalysis)
                .GetMethod("GetCutImpact", 
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            // Act - Test with the first cut
            var impact = (double)methodInfo.Invoke(analysis, new object[] { 0 });

            // Assert - Impact should be a non-negative value
            Assert.True(impact >= 0);
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
