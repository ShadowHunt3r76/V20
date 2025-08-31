using System;
using System.IO;
using Xunit;
using LinearProgramming.Algorithms.Knapsack;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Parsing;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Tests.SensitivityAnalysis
{
    public class KnapsackSensitivityAnalysisTests : IDisposable
    {
        private readonly KnapsackSolution _testSolution;
        private readonly ParsedLinearProgrammingModel _testModel;

        public KnapsackSensitivityAnalysisTests()
        {
            // Create a test model (0-1 Knapsack problem)
            _testModel = new ParsedLinearProgrammingModel
            {
                Objective = new Objective { Type = OptimizationType.Maximize },
                Variables = new[]
                {
                    new Variable { Name = "x1", Type = VariableType.Binary },
                    new Variable { Name = "x2", Type = VariableType.Binary },
                    new Variable { Name = "x3", Type = VariableType.Binary }
                },
                ObjectiveCoefficients = new[] { 10.0, 20.0, 15.0 },
                ConstraintCoefficients = new[] { new[] { 5.0, 10.0, 7.0 } },
                ConstraintConstants = new[] { 12.0 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual }
            };

            // Create a test solution
            _testSolution = new KnapsackSolution
            {
                BestObjectiveValue = 25.0,  // x1=1, x2=0, x3=1 (value=25, weight=12)
                BestSolution = new[] { 1.0, 0.0, 1.0 },
                OverallStatus = SolutionStatus.Optimal,
                TotalItems = 3,
                KnapsackCapacity = 12,
                Items = new[]
                {
                    new KnapsackItem { Value = 10, Weight = 5, IsSelected = true },
                    new KnapsackItem { Value = 20, Weight = 10, IsSelected = false },
                    new KnapsackItem { Value = 15, Weight = 7, IsSelected = true }
                },
                Model = _testModel
            };
        }

        [Fact]
        public void Constructor_WithValidParameters_InitializesCorrectly()
        {
            // Act
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);

            // Assert
            Assert.NotNull(analysis);
        }

        [Fact]
        public void PerformAnalysis_WithValidSolution_OutputsSensitivityInfo()
        {
            // Arrange
            var consoleOutput = new StringWriter();
            Console.SetOut(consoleOutput);
            
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);

            // Act
            analysis.PerformAnalysis();
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("KNAPSACK SENSITIVITY ANALYSIS", output);
            Assert.Contains("SOLUTION INFORMATION", output);
            Assert.Contains("ITEM SENSITIVITY", output);
            Assert.Contains("CAPACITY SENSITIVITY", output);
        }

        [Fact]
        public void PerformItemSensitivityAnalysis_WithValidSolution_OutputsItemInfo()
        {
            // Arrange
            var consoleOutput = new StringWriter();
            Console.SetOut(consoleOutput);
            
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);
            
            // Use reflection to access the protected method for testing
            var methodInfo = typeof(KnapsackSensitivityAnalysis)
                .GetMethod("PerformItemSensitivityAnalysis", 
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            // Act
            methodInfo.Invoke(analysis, null);
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("ITEM SENSITIVITY", output);
            Assert.Contains("Item", output);
            Assert.Contains("Value", output);
            Assert.Contains("Weight", output);
        }

        [Fact]
        public void GetValueSensitivity_WithValidItem_ReturnsCorrectSensitivity()
        {
            // Arrange
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);
            
            // Use reflection to access the protected method for testing
            var methodInfo = typeof(KnapsackSensitivityAnalysis)
                .GetMethod("GetValueSensitivity", 
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            // Act - Test with an included item (x1)
            var sensitivityIncluded = (double)methodInfo.Invoke(analysis, new object[] { 0 });
            
            // Test with an excluded item (x2)
            var sensitivityExcluded = (double)methodInfo.Invoke(analysis, new object[] { 1 });

            // Assert
            // For included items, sensitivity is how much the value can decrease before being excluded
            Assert.True(sensitivityIncluded > 0);
            // For excluded items, sensitivity is how much the value must increase to be included
            Assert.True(sensitivityExcluded > 0);
        }

        [Fact]
        public void GetCapacitySensitivity_WithValidSolution_ReturnsValidRange()
        {
            // Arrange
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);
            
            // Use reflection to access the protected method for testing
            var methodInfo = typeof(KnapsackSensitivityAnalysis)
                .GetMethod("GetCapacitySensitivity", 
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            // Act
            var range = ((double min, double max))methodInfo.Invoke(analysis, null);

            // Assert - The range should be reasonable values around the current capacity (12)
            Assert.True(range.min < 12.0);
            Assert.True(range.max > 12.0);
        }

        [Fact]
        public void GetObjectiveCoefficientRange_WithValidVariable_ReturnsValidRange()
        {
            // Arrange
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);

            // Act - Test with x1 (included in solution)
            var rangeIncluded = analysis.GetObjectiveCoefficientRange(0);
            
            // Test with x2 (excluded from solution)
            var rangeExcluded = analysis.GetObjectiveCoefficientRange(1);

            // Assert - Included item's value can decrease by its sensitivity
            Assert.True(rangeIncluded.lowerBound < 10.0);
            Assert.True(rangeIncluded.upperBound >= 10.0);
            
            // Excluded item's value needs to increase by its sensitivity to be included
            Assert.True(rangeExcluded.lowerBound <= 20.0);
            Assert.True(rangeExcluded.upperBound > 20.0);
        }

        [Fact]
        public void GetRightHandSideRange_WithValidConstraint_ReturnsValidRange()
        {
            // Arrange
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);

            // Act - Test with the capacity constraint
            var range = analysis.GetRightHandSideRange(0);

            // Assert - The range should be reasonable values around the current capacity (12)
            Assert.True(range.lowerBound < 12.0);
            Assert.True(range.upperBound > 12.0);
        }

        [Fact]
        public void GetDualValues_WithValidSolution_ReturnsCorrectDuals()
        {
            // Arrange
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var duals = analysis.GetDualValues();

            // Assert - Should return dual values for the capacity constraint
            Assert.Single(duals); // Only one constraint in knapsack
            // Dual value should be positive (capacity constraint is binding)
            Assert.True(duals[0] > 0);
        }

        [Fact]
        public void GetReducedCosts_WithValidSolution_ReturnsCorrectReducedCosts()
        {
            // Arrange
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var reducedCosts = analysis.GetReducedCosts();

            // Assert - Should return reduced costs for all items
            Assert.Equal(3, reducedCosts.Length); // Three items in our test case
            
            // For included items, reduced cost is negative of the value sensitivity
            Assert.True(reducedCosts[0] <= 0); // x1 is included
            // For excluded items, reduced cost is the value sensitivity
            Assert.True(reducedCosts[1] >= 0); // x2 is excluded
            Assert.True(reducedCosts[2] <= 0); // x3 is included
        }

        [Fact]
        public void IsOptimal_WithValidSolution_ReturnsTrue()
        {
            // Arrange
            var analysis = new KnapsackSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var isOptimal = analysis.IsOptimal();

            // Assert
            Assert.True(isOptimal);
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
