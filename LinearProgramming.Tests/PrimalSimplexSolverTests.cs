using System;
using Xunit;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Parsing;
using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;

namespace LinearProgramming.Tests
{
    public class PrimalSimplexSolverTests
    {
        private readonly PrimalSimplexSolver _solver;
        private readonly ParsedLinearProgrammingModel.UniversalLinearProgrammingParser _parser;

        public PrimalSimplexSolverTests()
        {
            _solver = new PrimalSimplexSolver();
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
                    new double[] { 1, 1, 1 },
                    new double[] { 2, 1, 0 },
                    new double[] { 0, 1, 1 }
                },
                RHSVector = new double[] { 4, 5, 3 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual, ConstraintType.LessThanOrEqual },
                ObjectiveCoefficients = new double[] { 3, 2, 4 },
                VariableTypes = new[] { VariableType.NonNegative, VariableType.NonNegative, VariableType.NonNegative },
                VariableNames = new[] { "x1", "x2", "x3" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue > 0);
            Assert.Equal(3, solution.VariableValues.Count);
        }

        [Fact]
        public void Solve_WithInfeasibleProblem_ReturnsInfeasibleStatus()
        {
            // Arrange - Infeasible constraints: x1 + x2 <= 1 and x1 + x2 >= 2
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, 1 },
                    new double[] { 1, 1 }
                },
                RightHandSide = new double[] { 1, 2 },
                ConstraintTypes = new[] { ConstraintType.LessThanOrEqual, ConstraintType.GreaterThanOrEqual },
                ObjectiveCoefficients = new double[] { 1, 1 },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Infeasible", solution.Status);
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
                ObjectiveCoefficients = new double[] { 1, 1 },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Unbounded", solution.Status);
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
                ObjectiveCoefficients = new double[] { 1, 1 },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue > 0);
        }

        [Fact]
        public void Solve_WithNegativeRHS_HandlesCorrectly()
        {
            // Arrange - Problem with negative RHS
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, 1 },
                    new double[] { -1, -1 }
                },
                RightHandSide = new double[] { -1, -2 },
                ConstraintTypes = new[] { ConstraintType.GreaterThanOrEqual, ConstraintType.GreaterThanOrEqual },
                ObjectiveCoefficients = new double[] { 1, 1 },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue >= 0);
        }

        [Fact]
        public void Solve_WithEqualityConstraints_HandlesCorrectly()
        {
            // Arrange - Problem with equality constraint
            var model = new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = new[]
                {
                    new double[] { 1, 1 },
                    new double[] { 1, -1 }
                },
                RightHandSide = new double[] { 1, 0 },
                ConstraintTypes = new[] { ConstraintType.Equal, ConstraintType.LessThanOrEqual },
                ObjectiveCoefficients = new double[] { 1, 1 },
                VariableNames = new[] { "x1", "x2" }
            };

            // Act
            var solution = _solver.Solve(model);

            // Assert
            Assert.Equal("Optimal", solution.Status);
            Assert.True(solution.ObjectiveValue > 0);
        }
    }
}
