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

            // Create variable names (x1, x2, ..., s1, s2, ..., a1, a2, ...)
            var varNames = new List<string>();
            for (int i = 0; i < n; i++)
                varNames.Add($"x{i + 1}");
            for (int i = 0; i < slackCount; i++)
                varNames.Add($"s{i + 1}");
            for (int i = 0; i < artificialCount; i++)
                varNames.Add($"a{i + 1}");

            // Get constraint types for the model
            var constraintTypes = Constraints.Select(con => con.Type).ToArray();

            return new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = A.Select(r => r.ToArray()).ToArray(),
                RHSVector = b.ToArray(),
                ObjectiveCoefficients = c.ToArray(),
                VariableTypes = canonicalVarTypes.ToArray(),
                ConstraintTypes = constraintTypes,
                VariableNames = varNames.ToArray()
            };
        }

        public class CanonicalLinearProgrammingModel
        {
            /// <summary>
            /// Gets or sets the coefficient matrix (A) for the constraints
            /// </summary>
            public double[][] CoefficientMatrix { get; set; }

            /// <summary>
            /// Gets or sets the right-hand side (RHS) values (b) for the constraints
            /// </summary>
            public double[] RHSVector { get; set; }

            /// <summary>
            /// Gets or sets the right-hand side values (b) for the constraints
            /// This is an alias for RHSVector for backward compatibility
            /// </summary>
            public double[] RightHandSide 
            { 
                get => RHSVector;
                set => RHSVector = value;
            }

            /// <summary>
            /// Gets or sets the objective function coefficients (c)
            /// </summary>
            public double[] ObjectiveCoefficients { get; set; }

            /// <summary>
            /// Gets or sets the types of variables in the model
            /// </summary>
            public VariableType[] VariableTypes { get; set; }

            /// <summary>
            /// Gets or sets the constraint types (≤, ≥, or =)
            /// </summary>
            public ConstraintType[] ConstraintTypes { get; set; }

            /// <summary>
            /// Gets or sets the names of the variables
            /// </summary>
            public string[] VariableNames { get; set; }
        }

        /// <summary>
        /// Universal parser for linear programming input files and strings.
        /// Supports both LP and IP models with the exact format specification.
        /// </summary>
        public class UniversalLinearProgrammingParser
        {
            // Performance optimization: Compile regex patterns once
            private static readonly Regex ObjectiveHeaderRegex = new Regex(@"^(max|min)(\s+[+-]\s*\d+(?:\.\d+)?)+$", 
                RegexOptions.Compiled | RegexOptions.IgnoreCase);
                
            private static readonly Regex ConstraintRegex = new Regex(
                @"^([+-]\s*\d+(?:\.\d+)?(?:\s+[+-]\s*\d+(?:\.\d+)?)*)\s*(<=|>=|=)\s*([+-]?\d+(?:\.\d+)?)\s*$", 
                RegexOptions.Compiled);
                
            private static readonly Regex VariableTypeRegex = new Regex(
                @"^(\+|-|urs|int|bin)(\s+(\+|-|urs|int|bin))*\s*$", 
                RegexOptions.Compiled | RegexOptions.IgnoreCase);
                
            private static readonly Regex WhitespaceRegex = new Regex("\\s+", 
                RegexOptions.Compiled);
                
            // Culture settings for number parsing
            private readonly IFormatProvider _numberFormat;
            private readonly NumberStyles _numberStyle;
            
            // Error messages
            private const string InvalidObjectiveFormat = "Invalid objective format. Expected: 'max/min [+-]coeff1 [+-]coeff2 ...' (e.g., 'max +2 +3 -1')";
            private const string InvalidConstraintFormat = "Invalid constraint format. Expected: '[+-]coeff1 [+-]coeff2 ... <=/>=/= rhs' (e.g., '+1 -2 +3 <= 4')";
            private const string InvalidVariableTypeFormat = "Invalid variable type specification. Expected: '+ - urs int bin' (one per variable)";
            private const string MismatchedVariableCount = "Number of variables in objective, constraints, and types must match";
            
            /// <summary>
            /// Initializes a new instance of the UniversalLinearProgrammingParser class.
            /// </summary>
            /// <param name="useInvariantCulture">Set to true to use invariant culture for number parsing (recommended for file I/O)</param>
            public UniversalLinearProgrammingParser(bool useInvariantCulture = true)
            {
                _numberFormat = useInvariantCulture ? CultureInfo.InvariantCulture : CultureInfo.CurrentCulture;
                _numberStyle = NumberStyles.Float | NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint;
            }
            // Regex patterns and error messages are already defined in the class

            /// <summary>
            /// Parses a linear programming model from a string.
            /// </summary>
            /// <param name="input">Input string containing the model</param>
            /// <returns>Parsed linear programming model</returns>
            /// <exception cref="FormatException">If the input format is invalid</exception>
            public ParsedLinearProgrammingModel ParseFromString(string input)
            {
                if (string.IsNullOrWhiteSpace(input))
                    throw new FormatException("Input cannot be empty");

                try
                {
                    string[] lines = input.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries)
                        .Select(line => line.Trim())
                        .Where(line => !string.IsNullOrWhiteSpace(line))
                        .ToArray();

                    if (lines.Length < 3) // At least objective, one constraint, and variable types
                        throw new FormatException("Input must contain at least 3 lines (objective, constraints, variable types)");

                    return ParseLines(lines);
                }
                catch (Exception ex) when (!(ex is FormatException))
                {
                    throw new FormatException($"Error parsing input: {ex.Message}", ex);
                }
            }

            /// <summary>
            /// Parses the entire model from an array of lines.
            /// </summary>
            private ParsedLinearProgrammingModel ParseLines(string[] lines)
            {
                var model = new ParsedLinearProgrammingModel
                {
                    Objective = ParseObjectiveLine(lines[0]),
                    Constraints = new List<ParsedConstraint>()
                };

                // Parse constraints (all lines except first and last)
                for (int i = 1; i < lines.Length - 1; i++)
                {
                    var constraint = ParseConstraintLine(lines[i]);
                    if (constraint.Coefficients.Count != model.Objective.Coefficients.Count)
                        throw new FormatException($"Constraint {i} has {constraint.Coefficients.Count} coefficients, expected {model.Objective.Coefficients.Count}");
                    model.Constraints.Add(constraint);
                }

                // Parse variable types (last line)
                ParseVariableTypes(lines[^1], model);

                // Verify variable counts match
                if (model.Objective.Coefficients.Count != model.Variables.Count)
                    throw new FormatException(MismatchedVariableCount);

                return model;
            }

            /// <summary>
            /// Parses a constraint line.
            /// Format: "[+-]coeff1 [+-]coeff2 ... <=/>=/= rhs"
            /// </summary>
            private ParsedConstraint ParseConstraintLine(string line)
            {
                if (string.IsNullOrWhiteSpace(line))
                    throw new FormatException(InvalidConstraintFormat);

                // Performance: Use span for efficient string manipulation
                ReadOnlySpan<char> lineSpan = line.Trim();
                
                // Find the operator (<=, >=, or =)
                int opPos = -1;
                string op = null;
                
                // Check for each operator type
                int lePos = lineSpan.IndexOf("<=".AsSpan());
                int gePos = lineSpan.IndexOf(">=".AsSpan());
                int eqPos = lineSpan.IndexOf('=');
                
                if (lePos >= 0) op = "<=";
                else if (gePos >= 0) op = ">=";
                else if (eqPos >= 0) op = "=";
                
                if (op == null)
                    throw new FormatException($"No valid constraint operator (<=, >=, =) found in: {line}");
                
                opPos = lineSpan.IndexOf(op, StringComparison.Ordinal);
                
                // Split into coefficients and RHS
                var coeffsPart = lineSpan.Slice(0, opPos).Trim();
                var rhsPart = lineSpan.Slice(opPos + op.Length).Trim();
                
                // Parse coefficients using span for better performance
                var coeffs = new List<double>();
                foreach (var coeffStr in coeffsPart.ToString().Split(' ', StringSplitOptions.RemoveEmptyEntries))
                {
                    if (!double.TryParse(coeffStr, _numberStyle, _numberFormat, out double coeff))
                        throw new FormatException($"Invalid coefficient format: {coeffStr}");
                    coeffs.Add(coeff);
                }

                // Parse RHS
                if (!double.TryParse(rhsPart, _numberStyle, _numberFormat, out double rhs))
                    throw new FormatException($"Invalid RHS value: {rhsPart}");

                // Parse constraint type
                var constraintType = op switch
                {
                    "<=" => ConstraintType.LessThanOrEqual,
                    ">=" => ConstraintType.GreaterThanOrEqual,
                    "=" => ConstraintType.Equal,
                    _ => throw new FormatException($"Invalid constraint operator: {op}")
                };

                return new ParsedConstraint
                {
                    Coefficients = coeffs,
                    Type = constraintType,
                    RHS = rhs
                };
            }

            /// <summary>
            /// Parses the variable types line (last line of input).
            /// Format: "+ - urs int bin" (one per variable)
            /// </summary>
            private void ParseVariableTypes(string line, ParsedLinearProgrammingModel model)
            {
                if (string.IsNullOrWhiteSpace(line))
                    throw new FormatException(InvalidVariableTypeFormat);

                // Performance: Use span for efficient string manipulation
                var typeSpans = line.Trim().AsSpan();
                var variables = new List<ParsedVariable>();
                int varIndex = 0;
                
                // Process each character in the line
                int start = 0;
                for (int i = 0; i <= typeSpans.Length; i++)
                {
                    if (i == typeSpans.Length || char.IsWhiteSpace(typeSpans[i]))
                    {
                        if (i > start)
                        {
                            var typeSpan = typeSpans.Slice(start, i - start);
                            var variable = new ParsedVariable
                            {
                                Index = varIndex++,
                                Name = $"x{varIndex}",
                                Type = typeSpan switch
                                {
                                    _ when typeSpan.Equals("+".AsSpan(), StringComparison.Ordinal) => VariableType.NonNegative,
                                    _ when typeSpan.Equals("-".AsSpan(), StringComparison.Ordinal) => VariableType.NonPositive,
                                    _ when typeSpan.Equals("urs".AsSpan(), StringComparison.OrdinalIgnoreCase) => VariableType.Unrestricted,
                                    _ when typeSpan.Equals("int".AsSpan(), StringComparison.OrdinalIgnoreCase) => VariableType.Integer,
                                    _ when typeSpan.Equals("bin".AsSpan(), StringComparison.OrdinalIgnoreCase) => VariableType.Binary,
                                    _ => throw new FormatException($"Invalid variable type: {typeSpan.ToString()}. Must be one of: +, -, urs, int, bin")
                                }
                            };
                            variables.Add(variable);
                        }
                        start = i + 1;
                    }
                }

                if (variables.Count == 0)
                    throw new FormatException(InvalidVariableTypeFormat);

                model.Variables = variables;
            }

            /// <summary>
            /// Parses the objective function line (first line of input).
            /// Format: "max/min [+-]coeff1 [+-]coeff2 ..."
            /// </summary>
            private ParsedObjectiveFunction ParseObjectiveLine(string line)
            {
                if (string.IsNullOrWhiteSpace(line))
                    throw new FormatException(InvalidObjectiveFormat);

                // Performance: Use span for efficient string manipulation
                ReadOnlySpan<char> lineSpan = line.Trim();
                
                // Find the first space to separate 'max/min' from coefficients
                int firstSpace = lineSpan.IndexOf(' ');
                if (firstSpace < 0)
                    throw new FormatException(InvalidObjectiveFormat);

                // Parse optimization type
                var optTypeSpan = lineSpan.Slice(0, firstSpace);
                if (!(optTypeSpan.Equals("max".AsSpan(), StringComparison.OrdinalIgnoreCase) || 
                      optTypeSpan.Equals("min".AsSpan(), StringComparison.OrdinalIgnoreCase)))
                {
                    throw new FormatException(InvalidObjectiveFormat);
                }
                var optType = optTypeSpan.Equals("max".AsSpan(), StringComparison.OrdinalIgnoreCase) 
                    ? OptimizationType.Maximize 
                    : OptimizationType.Minimize;

                // Parse coefficients using span for better performance
                var coeffSpan = lineSpan.Slice(firstSpace + 1);
                var coeffs = new List<double>();
                
                // Split by whitespace and parse each coefficient
                var coeffParts = coeffSpan.Trim().ToString().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                foreach (var part in coeffParts)
                {
                    if (!double.TryParse(part, _numberStyle, _numberFormat, out double coeff))
                        throw new FormatException($"Invalid coefficient format: {part}");
                    coeffs.Add(coeff);
                }

                if (coeffs.Count == 0)
                    throw new FormatException("At least one coefficient is required");

                return new ParsedObjectiveFunction
                {
                    Optimization = optType,
                    Coefficients = coeffs
                };
            }

            /// <summary>
            /// Parses a linear programming model from a file.
            /// </summary>
            /// <param name="filePath">Path to the input file</param>
            /// <returns>Parsed linear programming model</returns>
            /// <exception cref="FileNotFoundException">If the file does not exist</exception>
            /// <exception cref="FormatException">If the file format is invalid</exception>
            public ParsedLinearProgrammingModel ParseFromFile(string filePath)
            {
                if (!File.Exists(filePath))
                    throw new FileNotFoundException($"Input file not found: {filePath}");

                try
                {
                    string[] lines = File.ReadAllLines(filePath)
                        .Where(line => !string.IsNullOrWhiteSpace(line))
                        .Select(line => line.Trim())
                        .ToArray();

                    if (lines.Length < 3) // At least objective, one constraint, and variable types
                        throw new FormatException("Input file must contain at least 3 lines (objective, constraints, variable types)");

                    return ParseLines(lines);
                }
                catch (Exception ex) when (!(ex is FormatException))
                {
                    throw new FormatException($"Error parsing file: {ex.Message}", ex);
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