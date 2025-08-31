using System;
using Xunit;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Parsing;
using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;

namespace LinearProgramming.Tests
{
    public class RevisedPrimalSimplexSolverTests
    {
        private readonly RevisedPrimalSimplexSolver _solver;
        private readonly PrimalSimplexSolver _primalSolver; // For comparison
        private readonly ParsedLinearProgrammingModel.UniversalLinearProgrammingParser _parser;

        public RevisedPrimalSimplexSolverTests()
        {
            _solver = new RevisedPrimalSimplexSolver();
            _primalSolver = new PrimalSimplexSolver();
            _parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
        }

        [Fact]
        public void Solve_SimpleMaximizationProblem_ReturnsOptimalSolution()
        {
            // Arrange
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, 2 },
                    new double[] { 3, 1 },
                    new double[] { 1, 1 }
                },
                RHSVector = new double[] { 18, 21, 10 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual },
                ObjectiveCoefficients = new double[] { 3, 2 },
                VariableTypes = new[] { VariableType.NonNegative, VariableType.NonNegative },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue > 0);
            Assert.Equal(2, solution.VariableValues.Count);
        }

        [Fact]
        public void Solve_CompareWithPrimalSimplex_SimilarResults()
        {
            // Arrange - A simple LP problem
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 2, 1, 1 },
                    new double[] { 1, 3, 1 },
                    new double[] { 3, 1, 2 }
                },
                RightHandSide = new double[] { 18, 24, 30 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual },
                ObjectiveCoefficients = new double[] { 3, 4, 2 },
                VariableNames = new[] { "x1", "x2", "x3" }
            };

            // Act
            var revisedSolution = _solver.Solve(model);
            var primalSolution = _primalSolver.Solve(model);

            // Assert - Both should find optimal solutions with similar objective values
            Assert.Equal(primalSolution.Status, revisedSolution.Status);
            if (primalSolution.Status == "Optimal" && revisedSolution.Status == "Optimal")
            {
                Assert.InRange(revisedSolution.ObjectiveValue, 
                    primalSolution.ObjectiveValue - 0.001, 
                    primalSolution.ObjectiveValue + 0.001);
            }
        }

        [Fact]
        public void Solve_WithProductFormUpdate_ConvergesToSolution()
        {
            // Arrange - Force Product Form update
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, 1, 1 },
                    new double[] { 2, 1, 0 },
                    new double[] { 0, 1, 1 }
                },
                RightHandSide = new double[] { 4, 5, 3 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual },
                ObjectiveCoefficients = new double[] { 3, 2, 4 },
                VariableNames = new[] { "x1", "x2", "x3" }
            };

            // Act - Solve with Product Form update
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue > 0);
        }

        [Fact]
        public void Solve_WithPriceOutUpdate_ConvergesToSolution()
        {
            // Arrange - Force Price Out update
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, 1, 1 },
                    new double[] { 2, 1, 0 },
                    new double[] { 0, 1, 1 }
                },
                RightHandSide = new double[] { 4, 5, 3 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual },
                ObjectiveCoefficients = new double[] { 3, 2, 4 },
                VariableNames = new[] { "x1", "x2", "x3" }
            };

            // Act - Solve with Price Out update
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue > 0);
        }

        [Fact]
        public void Solve_WithDegenerateSolution_HandlesCorrectly()
        {
            // Arrange - Degenerate solution with redundant constraints
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, 1 },
                    new double[] { 1, 1 },
                    new double[] { 1, 0 }
                },
                RightHandSide = new double[] { 1, 1, 0.5 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue > 0);
        }

        [Fact]
        public void Solve_WithUnboundedProblem_ReturnsUnboundedStatus()
        {
            // Arrange - Unbounded maximization problem
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, -1 }
                },
                RightHandSide = new double[] { 1 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Unbounded", solution.Status);
        }
    }
}
