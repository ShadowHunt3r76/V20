using LinearProgramming.Algorithms;
using LinearProgramming.Parsing;
using System.Text.RegularExpressions;

namespace LinearProgramming.Tests;

public class ParserFormatTests
{
    [Fact]
    public void TestExactFormatFromSpecification()
    {
        // Test the exact format from the specification:
        // max +2 +3 +3 +5 +2 +4
        // +11 +8 +6 +14 +10 +10 <= 40
        // bin bin bin bin bin bin
        
        var parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
        string input = @"max +2 +3 +3 +5 +2 +4
+11 +8 +6 +14 +10 +10 <= 40
bin bin bin bin bin bin";
        
        // This should not throw an exception
        var parsedModel = parser.ParseFromString(input);
        
        // Verify the parsing worked correctly
        Assert.NotNull(parsedModel);
        Assert.Equal(OptimizationType.Maximize, parsedModel.Objective.Optimization);
        Assert.Equal(6, parsedModel.Objective.Coefficients.Count);
        Assert.Equal(2.0, parsedModel.Objective.Coefficients[0]);
        Assert.Equal(3.0, parsedModel.Objective.Coefficients[1]);
        Assert.Equal(4.0, parsedModel.Objective.Coefficients[5]);
        
        Assert.Single(parsedModel.Constraints);
        Assert.Equal(ConstraintType.LessThanOrEqual, parsedModel.Constraints[0].Type);
        Assert.Equal(40.0, parsedModel.Constraints[0].RHS);
        Assert.Equal(6, parsedModel.Constraints[0].Coefficients.Count);
        Assert.Equal(11.0, parsedModel.Constraints[0].Coefficients[0]);
        Assert.Equal(10.0, parsedModel.Constraints[0].Coefficients[5]);
        
        Assert.Equal(6, parsedModel.Variables.Count);
        Assert.All(parsedModel.Variables, v => Assert.Equal(VariableType.Binary, v.Type));
    }
    
    [Fact]
    public void TestMultipleConstraintsFormat()
    {
        // Test with multiple constraints as mentioned in requirements
        var parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
        string input = @"max +2 +3 +1
+1 +2 +0 <= 4
+2 +1 +3 >= 5
+1 +1 +1 = 3
+ + +";
        
        var parsedModel = parser.ParseFromString(input);
        
        Assert.NotNull(parsedModel);
        Assert.Equal(OptimizationType.Maximize, parsedModel.Objective.Optimization);
        Assert.Equal(3, parsedModel.Objective.Coefficients.Count);
        
        Assert.Equal(3, parsedModel.Constraints.Count);
        Assert.Equal(ConstraintType.LessThanOrEqual, parsedModel.Constraints[0].Type);
        Assert.Equal(ConstraintType.GreaterThanOrEqual, parsedModel.Constraints[1].Type);
        Assert.Equal(ConstraintType.Equal, parsedModel.Constraints[2].Type);
        
        Assert.Equal(3, parsedModel.Variables.Count);
        Assert.All(parsedModel.Variables, v => Assert.Equal(VariableType.NonNegative, v.Type));
    }
    
    [Fact]
    public void TestConstraintRegexPattern()
    {
        // Test the constraint regex pattern directly
        var constraintRegex = new Regex(@"^([+-]?\d+(?:\.\d+)?(?:\s+[+-]?\d+(?:\.\d+)?)*\s*)(<=|>=|=)\s*([+-]?\d+(?:\.\d+)?)$", RegexOptions.IgnoreCase);
        
        string testConstraint = "+11 +8 +6 +14 +10 +10 <= 40";
        var match = constraintRegex.Match(testConstraint);
        
        Assert.True(match.Success, $"Constraint regex should match: '{testConstraint}'");
        Assert.Equal("+11 +8 +6 +14 +10 +10 ", match.Groups[1].Value);
        Assert.Equal("<=", match.Groups[2].Value);
        Assert.Equal("40", match.Groups[3].Value);
    }
    
    [Fact]
    public void TestCoefficientParsing()
    {
        // Test coefficient parsing with both positive and negative values
        var parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
        string input = @"min -2 +3 -1
+1 -2 +0 <= 4
+ - +";
        
        var parsedModel = parser.ParseFromString(input);
        
        Assert.Equal(OptimizationType.Minimize, parsedModel.Objective.Optimization);
        Assert.Equal(-2.0, parsedModel.Objective.Coefficients[0]);
        Assert.Equal(3.0, parsedModel.Objective.Coefficients[1]);
        Assert.Equal(-1.0, parsedModel.Objective.Coefficients[2]);
        
        Assert.Equal(1.0, parsedModel.Constraints[0].Coefficients[0]);
        Assert.Equal(-2.0, parsedModel.Constraints[0].Coefficients[1]);
        Assert.Equal(0.0, parsedModel.Constraints[0].Coefficients[2]);
        
        Assert.Equal(VariableType.NonNegative, parsedModel.Variables[0].Type);
        Assert.Equal(VariableType.NonPositive, parsedModel.Variables[1].Type);
        Assert.Equal(VariableType.NonNegative, parsedModel.Variables[2].Type);
    }
}
