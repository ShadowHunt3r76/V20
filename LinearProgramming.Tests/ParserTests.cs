using System;
using System.IO;
using Xunit;
using LinearProgramming.Parsing;

namespace LinearProgramming.Tests
{
    public class ParserTests
    {
        private readonly ParsedLinearProgrammingModel.UniversalLinearProgrammingParser _parser;

        public ParserTests()
        {
            _parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
        }

        [Fact]
        public void ParseFromString_ValidMaximizationProblem_ReturnsCorrectModel()
        {
            // Arrange
            string input = "max +2 +3 +4\n" +
                         "+1 +2 +3 <= 10\n" +
                         "+4 +5 +6 >= 20\n" +
                         "+7 +8 +9 = 30\n" +
                         "+ + +";

            // Act
            var model = _parser.ParseFromString(input);

            // Assert
            Assert.Equal(OptimizationType.Maximize, model.Objective.Optimization);
            Assert.Equal(new[] { 2.0, 3.0, 4.0 }, model.Objective.Coefficients);
            Assert.Equal(3, model.Constraints.Count);
            Assert.Equal(new[] { 1.0, 2.0, 3.0 }, model.Constraints[0].Coefficients);
            Assert.Equal(ConstraintType.LessThanOrEqual, model.Constraints[0].Type);
            Assert.Equal(10.0, model.Constraints[0].RHS);
            Assert.Equal(3, model.Variables.Count);
            Assert.All(model.Variables, v => Assert.Equal(VariableType.NonNegative, v.Type));
        }

        [Fact]
        public void ParseFromString_ValidMinimizationProblem_ReturnsCorrectModel()
        {
            // Arrange
            string input = "min -1.5 +2.25 -3.75\n" +
                         "-1 -2 -3 >= -10\n" +
                         "+1 -1 +1 <= 5\n" +
                         "+ + +";

            // Act
            var model = _parser.ParseFromString(input);

            // Assert
            Assert.Equal(OptimizationType.Minimize, model.Objective.Optimization);
            Assert.Equal(new[] { -1.5, 2.25, -3.75 }, model.Objective.Coefficients);
            Assert.Equal(2, model.Constraints.Count);
            Assert.Equal(new[] { -1.0, -2.0, -3.0 }, model.Constraints[0].Coefficients);
            Assert.Equal(ConstraintType.GreaterThanOrEqual, model.Constraints[0].Type);
            Assert.Equal(-10.0, model.Constraints[0].RHS);
        }

        [Fact]
        public void ParseFromString_WithDifferentVariableTypes_ReturnsCorrectModel()
        {
            // Arrange
            string input = "max +1 +2 +3 +4 +5\n" +
                         "+1 +2 +3 +4 +5 <= 100\n" +
                         "+ - urs int bin";

            // Act
            var model = _parser.ParseFromString(input);

            // Assert
            Assert.Equal(5, model.Variables.Count);
            Assert.Equal(VariableType.NonNegative, model.Variables[0].Type);
            Assert.Equal(VariableType.NonPositive, model.Variables[1].Type);
            Assert.Equal(VariableType.Unrestricted, model.Variables[2].Type);
            Assert.Equal(VariableType.Integer, model.Variables[3].Type);
            Assert.Equal(VariableType.Binary, model.Variables[4].Type);
        }

        [Fact]
        public void ParseFromString_WithDecimalValues_ReturnsCorrectModel()
        {
            // Arrange
            string input = "max +1.5 -2.25 +3.75\n" +
                         "+0.1 -0.2 +0.3 <= 1.5\n" +
                         "+ + +";

            // Act
            var model = _parser.ParseFromString(input);

            // Assert
            Assert.Equal(new[] { 1.5, -2.25, 3.75 }, model.Objective.Coefficients);
            Assert.Equal(new[] { 0.1, -0.2, 0.3 }, model.Constraints[0].Coefficients);
            Assert.Equal(1.5, model.Constraints[0].RHS);
        }

        [Fact]
        public void ParseFromString_WithMultipleSpaces_ReturnsCorrectModel()
        {
            // Arrange
            string input = "max   +1   +2   +3\n  +4   +5   +6  <=  100  \n  +   -   urs";

            // Act
            var model = _parser.ParseFromString(input);

            // Assert
            Assert.Equal(new[] { 1.0, 2.0, 3.0 }, model.Objective.Coefficients);
            Assert.Equal(new[] { 4.0, 5.0, 6.0 }, model.Constraints[0].Coefficients);
            Assert.Equal(VariableType.NonNegative, model.Variables[0].Type);
            Assert.Equal(VariableType.NonPositive, model.Variables[1].Type);
            Assert.Equal(VariableType.Unrestricted, model.Variables[2].Type);
        }

        [Fact]
        public void ParseFromString_WithEmptyInput_ThrowsFormatException()
        {
            // Arrange
            string input = "";

            // Act & Assert
            var exception = Assert.Throws<FormatException>(() => _parser.ParseFromString(input));
            Assert.Contains("cannot be empty", exception.Message);
        }

        [Fact]
        public void ParseFromString_WithInvalidObjectiveFormat_ThrowsFormatException()
        {
            // Arrange
            string input = @"invalid +1 +2
+1 +2 <= 10
+ +";

            // Act & Assert
            var exception = Assert.Throws<FormatException>(() => _parser.ParseFromString(input));
            Assert.Contains("Invalid objective format", exception.Message);
        }

        [Fact]
        public void ParseFromString_WithMismatchedVariableCounts_ThrowsFormatException()
        {
            // Arrange - Objective has 3 vars, constraint has 2, types have 4
            string input = @"max +1 +2 +3
+1 +2 <= 10
+ + + +";

            // Act & Assert
            var exception = Assert.Throws<FormatException>(() => _parser.ParseFromString(input));
            Assert.Contains("Number of variables", exception.Message);
        }

        [Fact]
        public void ParseFromString_WithInvalidVariableType_ThrowsFormatException()
        {
            // Arrange
            string input = @"max +1 +2
+1 +2 <= 10
invalid invalid";

            // Act & Assert
            var exception = Assert.Throws<FormatException>(() => _parser.ParseFromString(input));
            Assert.Contains("Invalid variable type", exception.Message);
        }

        [Fact]
        public void ParseFromFile_ValidFile_ReturnsCorrectModel()
        {
            // Arrange
            string tempFile = Path.GetTempFileName();
            try
            {
                File.WriteAllText(tempFile, @"max +1 +2 +3
+1 +2 +3 <= 10
+4 +5 +6 >= 20
+7 +8 +9 = 30
+ + +");

                // Act
                var model = _parser.ParseFromFile(tempFile);

                // Assert
                Assert.Equal(OptimizationType.Maximize, model.Objective.Optimization);
                Assert.Equal(new[] { 1.0, 2.0, 3.0 }, model.Objective.Coefficients);
                Assert.Equal(3, model.Constraints.Count);
            }
            finally
            {
                if (File.Exists(tempFile))
                    File.Delete(tempFile);
            }
        }

        [Fact]
        public void ParseFromFile_NonExistentFile_ThrowsFileNotFoundException()
        {
            // Arrange
            string nonExistentFile = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString());

            // Act & Assert
            Assert.Throws<FileNotFoundException>(() => _parser.ParseFromFile(nonExistentFile));
        }

        [Fact]
        public void ToCanonicalForm_WithAllConstraintTypes_ReturnsCorrectCanonicalForm()
        {
            // Arrange
            var model = new ParsedLinearProgrammingModel
            {
                Objective = new ParsedObjectiveFunction
                {
                    Optimization = OptimizationType.Maximize,
                    Coefficients = new List<double> { 3, 2 }
                },
                Constraints = new List<ParsedConstraint>
                {
                    new ParsedConstraint { Coefficients = new List<double> { 2, 1 }, Type = ConstraintType.LessThanOrEqual, RHS = 100 },
                    new ParsedConstraint { Coefficients = new List<double> { 1, 1 }, Type = ConstraintType.GreaterThanOrEqual, RHS = 80 },
                    new ParsedConstraint { Coefficients = new List<double> { 1, 3 }, Type = ConstraintType.Equal, RHS = 90 }
                },
                Variables = new List<ParsedVariable>
                {
                    new ParsedVariable { Type = VariableType.NonNegative, Index = 0, Name = "x1" },
                    new ParsedVariable { Type = VariableType.NonNegative, Index = 1, Name = "x2" }
                }
            };

            // Act
            var canonical = model.ToCanonicalForm();

            // Assert
            // Should have 2 original variables + 1 slack + 1 surplus + 1 artificial for >= + 1 artificial for =
            Assert.Equal(5, canonical.ObjectiveCoefficients.Length);
            Assert.Equal(3, canonical.CoefficientMatrix.Length); // 3 constraints
            Assert.Equal(3, canonical.RHSVector.Length); // 3 constraints
            
            // Verify objective coefficients (original vars + slack/surplus/artificial vars)
            Assert.Equal(new[] { 3.0, 2.0, 0, 0, 0 }, canonical.ObjectiveCoefficients);
            
            // Verify variable types (original + slack/surplus/artificial)
            Assert.Equal(5, canonical.VariableTypes.Length);
            Assert.All(canonical.VariableTypes, t => 
                Assert.True(t == VariableType.NonNegative || t == VariableType.Integer || t == VariableType.Binary));
        }
    }
}
