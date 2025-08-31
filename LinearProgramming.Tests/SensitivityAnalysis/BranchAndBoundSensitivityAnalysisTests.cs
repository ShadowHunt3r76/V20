using System;
using System.IO;
using Xunit;
using LinearProgramming.Algorithms.BranchAndBound;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Parsing;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Tests.SensitivityAnalysis
{
    public class BranchAndBoundSensitivityAnalysisTests : IDisposable
    {
        private readonly BranchAndBoundSolution _testSolution;
        private readonly ParsedLinearProgrammingModel _testModel;

        public BranchAndBoundSensitivityAnalysisTests()
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
            _testSolution = new BranchAndBoundSolution
            {
                BestObjectiveValue = 30.0,  // x1=0, x2=1, x3=1 (value=35, but weight=17 > 12, so infeasible)
                BestSolution = new[] { 1.0, 1.0, 0.0 },  // x1=1, x2=1, x3=0 (value=30, weight=15 <= 12? Wait, 5+10=15 > 12, this is infeasible)
                // Let me correct this to a feasible solution
                BestObjectiveValue = 25.0,  // x1=1, x2=0, x3=1 (value=25, weight=12 <= 12)
                BestSolution = new[] { 1.0, 0.0, 1.0 },  // x1=1, x2=0, x3=1 (value=25, weight=12)
                OverallStatus = SolutionStatus.Optimal,
                TotalNodesExplored = 5,
                BestCandidate = new BranchAndBoundNode
                {
                    Solution = new PrimalSimplexSolution
                    {
                        Status = SolutionStatus.Optimal,
                        ObjectiveValue = 25.0,
                        VariableValues = new[] { 1.0, 0.0, 1.0 },
                        // Mock optimal tableau for the LP relaxation at the best node
                        OptimalTable = new double[,]
                        {
                            { 1.0, 0.0, 0.0, 2.0, 25.0 },  // Z + 2s1 = 25
                            { 0.0, 1.0, 0.0, 0.2, 1.0 },    // x1 + 0.2s1 = 1
                            { 0.0, 0.0, 1.0, -0.2, 1.0 }    // x3 - 0.2s1 = 1
                        }
                    },
                    Model = new ParsedLinearProgrammingModel
                    {
                        Variables = _testModel.Variables,
                        ObjectiveCoefficients = _testModel.ObjectiveCoefficients,
                        ConstraintTypes = _testModel.ConstraintTypes,
                        ConstraintCoefficients = _testModel.ConstraintCoefficients,
                        ConstraintConstants = _testModel.ConstraintConstants
                    }
                },
                CanonicalForm = new CanonicalForm
                {
                    ObjectiveCoefficients = _testModel.ObjectiveCoefficients,
                    ConstraintCoefficients = _testModel.ConstraintCoefficients,
                    ConstraintConstants = _testModel.ConstraintConstants,
                    ConstraintTypes = _testModel.ConstraintTypes
                }
            };
        }

        [Fact]
        public void Constructor_WithValidParameters_InitializesCorrectly()
        {
            // Act
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);

            // Assert
            Assert.NotNull(analysis);
        }

        [Fact]
        public void PerformAnalysis_WithValidSolution_OutputsSensitivityInfo()
        {
            // Arrange
            var consoleOutput = new StringWriter();
            Console.SetOut(consoleOutput);
            
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);

            // Act
            analysis.PerformAnalysis();
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("BRANCH AND BOUND SENSITIVITY ANALYSIS", output);
            Assert.Contains("SOLUTION INFORMATION", output);
            Assert.Contains("LP-BASED SENSITIVITY ANALYSIS", output);
            Assert.Contains("INTEGER-SPECIFIC SENSITIVITY ANALYSIS", output);
        }

        [Fact]
        public void PerformIntegerSensitivityAnalysis_WithValidSolution_OutputsIntegerSpecificInfo()
        {
            // Arrange
            var consoleOutput = new StringWriter();
            Console.SetOut(consoleOutput);
            
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);
            
            // Use reflection to access the protected method for testing
            var methodInfo = typeof(BranchAndBoundSensitivityAnalysis)
                .GetMethod("PerformIntegerSensitivityAnalysis", 
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            // Act
            methodInfo.Invoke(analysis, null);
            var output = consoleOutput.ToString();

            // Assert
            Assert.Contains("INTEGER-SPECIFIC SENSITIVITY ANALYSIS", output);
            Assert.Contains("Variable", output);
            Assert.Contains("Value", output);
        }

        [Fact]
        public void GetObjectiveCoefficientRange_WithValidVariable_ReturnsValidRange()
        {
            // Arrange
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);

            // Act - Test with a variable that's in the solution (x1)
            var range = analysis.GetObjectiveCoefficientRange(0);

            // Assert - The range should be reasonable values around the original coefficient (10.0)
            Assert.True(range.lowerBound < 10.0);
            Assert.True(range.upperBound > 10.0);
        }

        [Fact]
        public void GetRightHandSideRange_WithValidConstraint_ReturnsValidRange()
        {
            // Arrange
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);

            // Act - Test with the knapsack capacity constraint
            var range = analysis.GetRightHandSideRange(0);

            // Assert - The range should be reasonable values around the original RHS (12.0)
            Assert.True(range.lowerBound < 12.0);
            Assert.True(range.upperBound > 12.0);
        }

        [Fact]
        public void GetDualValues_WithValidSolution_ReturnsCorrectDuals()
        {
            // Arrange
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var duals = analysis.GetDualValues();

            // Assert - Should return dual values from the LP relaxation
            Assert.Single(duals); // Only one constraint in our test case
            Assert.True(duals[0] >= 0); // Dual should be non-negative for ≤ constraint
        }

        [Fact]
        public void GetReducedCosts_WithValidSolution_ReturnsCorrectReducedCosts()
        {
            // Arrange
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);

            // Act
            var reducedCosts = analysis.GetReducedCosts();

            // Assert - Should return reduced costs for all variables
            Assert.Equal(3, reducedCosts.Length); // Three variables in our test case
            // Basic variables should have reduced cost of 0 (x1 and x3 are in the basis)
            Assert.Equal(0.0, reducedCosts[0], 6); // x1 is basic
            Assert.Equal(0.0, reducedCosts[2], 6); // x3 is basic
        }

        [Fact]
        public void IsOptimal_WithValidSolution_ReturnsTrue()
        {
            // Arrange
            var analysis = new BranchAndBoundSensitivityAnalysis(_testSolution, _testModel);

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
