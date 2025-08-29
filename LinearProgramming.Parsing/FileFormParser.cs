using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace LinearProgramming.Parsing
{
    /// <summary>
    /// Specifies the optimization direction for the objective function.
    /// </summary>
    public enum OptimizationType
    {
        Maximize,
        Minimize
    }

    /// <summary>
    /// Specifies the relationship type for constraints.
    /// </summary>
    public enum ConstraintType
    {
        LessThanOrEqual,
        GreaterThanOrEqual,
        Equal
    }

    /// <summary>
    /// Specifies the type and bounds for variables.
    /// </summary>
    public enum VariableType
    {
        NonNegative, // +
        NonPositive, // -
        Unrestricted, // urs
        Integer,      // int
        Binary        // bin
    }

    /// <summary>
    /// Represents a parsed constraint with coefficients, type, and RHS value.
    /// </summary>
    public class ParsedConstraint
    {
        public List<double> Coefficients { get; set; }
        public ConstraintType Type { get; set; }
        public double RHS { get; set; }
    }

    /// <summary>
    /// Represents the objective function with optimization type and coefficients.
    /// </summary>
    public class ParsedObjectiveFunction
    {
        public OptimizationType Optimization { get; set; }
        public List<double> Coefficients { get; set; }
    }

    /// <summary>
    /// Represents a variable with type and metadata.
    /// </summary>
    public class ParsedVariable
    {
        public VariableType Type { get; set; }
        public int Index { get; set; }
        public string Name { get; set; }
    }

    /// <summary>
    /// Main model container for parsed linear programming data.
    /// </summary>
    public class ParsedLinearProgrammingModel
    {
        public ParsedObjectiveFunction Objective { get; set; }
        public List<ParsedConstraint> Constraints { get; set; }
        public List<ParsedVariable> Variables { get; set; }

        /// <summary>
        /// Returns the coefficient matrix (A) for constraints.
        /// </summary>
        public double[,] GetCoefficientMatrix()
        {
            int rows = Constraints.Count;
            int cols = Variables.Count;
            double[,] matrix = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    matrix[i, j] = Constraints[i].Coefficients[j];
            return matrix;
        }

        /// <summary>
        /// Returns the RHS vector (b) for constraints.
        /// </summary>
        public double[] GetRHSVector()
        {
            return Constraints.Select(c => c.RHS).ToArray();
        }

        /// <summary>
        /// Returns the objective coefficients (c).
        /// </summary>
        public double[] GetObjectiveCoefficients()
        {
            return Objective.Coefficients.ToArray();
        }

        /// <summary>
        /// Converts the parsed model to canonical (standard) form for simplex algorithms.
        /// Adds slack/surplus/artificial variables as needed.
        /// </summary>
        public CanonicalLinearProgrammingModel ToCanonicalForm()
        {
            // Determine number of original variables
            int n = Variables.Count;
            int m = Constraints.Count;
            var A = new List<List<double>>();
            var b = new List<double>();
            var c = new List<double>(Objective.Coefficients);
            var slackSurplus = new List<List<double>>();
            var canonicalVarTypes = new List<VariableType>(Variables.Select(v => v.Type));
            int slackCount = 0, artificialCount = 0;

            // Build canonical matrix
            for (int i = 0; i < m; i++)
            {
                var row = new List<double>(Constraints[i].Coefficients);
                var slackVars = new List<double>(new double[m]);
                var type = Constraints[i].Type;
                if (type == ConstraintType.LessThanOrEqual)
                {
                    slackVars[i] = 1;
                    slackCount++;
                    canonicalVarTypes.Add(VariableType.NonNegative);
                }
                else if (type == ConstraintType.GreaterThanOrEqual)
                {
                    slackVars[i] = -1;
                    slackCount++;
                    canonicalVarTypes.Add(VariableType.NonNegative);
                    // Add artificial variable for >=
                    row.AddRange(slackVars);
                    row.Add(1); // artificial
                    canonicalVarTypes.Add(VariableType.NonNegative);
                    artificialCount++;
                    A.Add(row);
                    b.Add(Constraints[i].RHS);
                    continue;
                }
                else if (type == ConstraintType.Equal)
                {
                    // Add artificial variable for =
                    slackVars[i] = 0;
                    row.AddRange(slackVars);
                    row.Add(1); // artificial
                    canonicalVarTypes.Add(VariableType.NonNegative);
                    artificialCount++;
                    A.Add(row);
                    b.Add(Constraints[i].RHS);
                    continue;
                }
                row.AddRange(slackVars);
                A.Add(row);
                b.Add(Constraints[i].RHS);
            }

            // Extend objective coefficients for slack/surplus/artificial variables
            int totalVars = n + m + artificialCount;
            while (c.Count < totalVars)
                c.Add(0);

            return new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = A.Select(r => r.ToArray()).ToArray(),
                RHSVector = b.ToArray(),
                ObjectiveCoefficients = c.ToArray(),
                VariableTypes = canonicalVarTypes.ToArray()
            };
        }

        public class CanonicalLinearProgrammingModel
        {
            public double[][] CoefficientMatrix { get; set; }
            public double[] RHSVector { get; set; }
            public double[] ObjectiveCoefficients { get; set; }
            public VariableType[] VariableTypes { get; set; }
        }

        /// <summary>
        /// Universal parser for linear programming input files and strings.
        /// Supports both LP and IP models.
        /// </summary>
        public class UniversalLinearProgrammingParser
        {
            // Updated regex patterns to handle the exact format specification
            private static readonly Regex ObjectiveRegex = new Regex(@"^(max|min)\s+([+-]\d+(?:\.\d+)?(?:\s+[+-]\d+(?:\.\d+)?)*)\s*$", RegexOptions.IgnoreCase);
            private static readonly Regex ConstraintRegex = new Regex(@"^([+-]\d+(?:\.\d+)?(?:\s+[+-]\d+(?:\.\d+)?)*)\s+(<=|>=|=)\s+([+-]?\d+(?:\.\d+)?)\s*$", RegexOptions.IgnoreCase);
            private static readonly Regex VariableTypeRegex = new Regex(@"^(\+|-|urs|int|bin)(\s+(\+|-|urs|int|bin))*$", RegexOptions.IgnoreCase);

            /// <summary>
            /// Parses a linear programming model from a text string.
            /// Throws detailed exceptions for format errors.
            /// </summary>
            public ParsedLinearProgrammingModel ParseFromString(string input)
            {
                if (string.IsNullOrWhiteSpace(input))
                    throw new LinearProgrammingParseException("Input string is empty.");

                var lines = input.Split(new string[] { "\n", "\r\n" }, StringSplitOptions.RemoveEmptyEntries)
                    .Select(l => l.Trim()).Where(l => !string.IsNullOrWhiteSpace(l)).ToList();
                if (lines.Count < 3)
                    throw new LinearProgrammingParseException("Input must have at least objective, one constraint, and variable types.");

                // Parse objective
                var objectiveLine = lines[0];
                var objectiveParts = objectiveLine.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                if (objectiveParts.Length < 2)
                    throw new LinearProgrammingParseException($"Objective function line is invalid: '{objectiveLine}'");
                if (objectiveParts[0].ToLower() != "max" && objectiveParts[0].ToLower() != "min")
                    throw new LinearProgrammingParseException($"Objective function must start with 'max' or 'min': '{objectiveLine}'");
                OptimizationType optType = objectiveParts[0].ToLower() == "max" ? OptimizationType.Maximize : OptimizationType.Minimize;
                var objCoeffs = ParseCoefficients(objectiveParts.Skip(1), 1);
                var objective = new ParsedObjectiveFunction { Optimization = optType, Coefficients = objCoeffs };

                // Parse constraints
                var constraints = new List<ParsedConstraint>();
                for (int i = 1; i < lines.Count - 1; i++)
                {
                    var constraintLine = lines[i];
                    
                    // Find the constraint operator (<=, >=, =)
                    string[] operators = { "<=", ">=", "=" };
                    string foundOperator = null;
                    int operatorIndex = -1;
                    
                    foreach (var op in operators)
                    {
                        int idx = constraintLine.IndexOf(op);
                        if (idx != -1)
                        {
                            foundOperator = op;
                            operatorIndex = idx;
                            break;
                        }
                    }
                    
                    if (foundOperator == null)
                        throw new LinearProgrammingParseException($"No valid constraint operator found on line {i + 1}: '{constraintLine}'");
                    
                    // Split the line into coefficients and RHS
                    string coeffPart = constraintLine.Substring(0, operatorIndex).Trim();
                    string rhsPart = constraintLine.Substring(operatorIndex + foundOperator.Length).Trim();
                    
                    // Parse coefficients
                    var coeffParts = coeffPart.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    var coeffs = ParseCoefficients(coeffParts, i + 1);
                    
                    // Parse constraint type
                    var type = foundOperator == "<=" ? ConstraintType.LessThanOrEqual :
                        foundOperator == ">=" ? ConstraintType.GreaterThanOrEqual : ConstraintType.Equal;
                    
                    // Parse RHS - accepts both formats like "40" and "+40"
                    if (!double.TryParse(rhsPart, NumberStyles.Float, CultureInfo.InvariantCulture, out double rhs))
                        throw new LinearProgrammingParseException($"RHS value format error on line {i + 1}: '{rhsPart}'");
                    
                    constraints.Add(new ParsedConstraint { Coefficients = coeffs, Type = type, RHS = rhs });
                }

                // Parse variable types
                var varTypeLine = lines.Last();
                var varTypeParts = varTypeLine.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                var variables = new List<ParsedVariable>();
                for (int i = 0; i < varTypeParts.Length; i++)
                {
                    try
                    {
                        variables.Add(new ParsedVariable
                        {
                            Type = ParseVariableType(varTypeParts[i]),
                            Index = i,
                            Name = $"x{i + 1}"
                        });
                    }
                    catch (FormatException ex)
                    {
                        throw new LinearProgrammingParseException($"Variable type format error at position {i + 1}: '{varTypeParts[i]}'", ex);
                    }
                }

                // Validation
                int varCount = objCoeffs.Count;
                if (constraints.Any(c => c.Coefficients.Count != varCount))
                    throw new LinearProgrammingParseException($"All constraints must have {varCount} coefficients matching the objective function.");
                if (variables.Count != varCount)
                    throw new LinearProgrammingParseException($"Number of variable types ({variables.Count}) must match number of objective coefficients ({varCount}).");

                return new ParsedLinearProgrammingModel
                {
                    Objective = objective,
                    Constraints = constraints,
                    Variables = variables
                };
            }

            /// <summary>
            /// Parses a linear programming model from a file.
            /// Throws detailed exceptions for format errors.
            /// </summary>
            public ParsedLinearProgrammingModel ParseFromFile(string filePath)
            {
                if (!File.Exists(filePath))
                    throw new LinearProgrammingParseException($"File not found: {filePath}");
                var content = File.ReadAllText(filePath);
                return ParseFromString(content);
            }

            private List<double> ParseCoefficients(IEnumerable<string> parts, int lineNumber)
            {
                var coeffs = new List<double>();
                int col = 1;
                foreach (var part in parts)
                {
                    if (string.IsNullOrWhiteSpace(part)) continue;
                    
                    // Accept coefficients with explicit signs (required for objective and constraint coefficients)
                    if (!Regex.IsMatch(part, @"^[+-]\d+(?:\.\d+)?$"))
                        throw new LinearProgrammingParseException($"Coefficient format error at line {lineNumber}, position {col}: '{part}' (coefficients must have explicit + or - signs)");
                    try
                    {
                        coeffs.Add(double.Parse(part, CultureInfo.InvariantCulture));
                    }
                    catch (Exception ex)
                    {
                        throw new LinearProgrammingParseException($"Coefficient parse error at line {lineNumber}, position {col}: '{part}'", ex);
                    }
                    col++;
                }
                return coeffs;
            }

            private VariableType ParseVariableType(string typeStr)
            {
                switch (typeStr.ToLower())
                {
                    case "+": return VariableType.NonNegative;
                    case "-": return VariableType.NonPositive;
                    case "urs": return VariableType.Unrestricted;
                    case "int": return VariableType.Integer;
                    case "bin": return VariableType.Binary;
                    default: throw new FormatException($"Unknown variable type: {typeStr}");
                }
            }
        }

        /// <summary>
        /// Custom exception for advanced error handling in LP parser.
        /// </summary>
        public class LinearProgrammingParseException : Exception
        {
            public LinearProgrammingParseException(string message) : base(message) { }
            public LinearProgrammingParseException(string message, Exception inner) : base(message, inner) { }
        }
    }
}