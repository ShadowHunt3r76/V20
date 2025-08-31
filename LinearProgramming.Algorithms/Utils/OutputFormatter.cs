using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using LinearProgramming.Parsing;

namespace LinearProgramming.Algorithms.Utils
{
    /// <summary>
    /// Provides consistent formatting for solver output across the application
    /// </summary>
    public static class OutputFormatter
    {
        // Box drawing characters
        private const char TopLeft = '╔';
        private const char TopRight = '╗';
        private const char BottomLeft = '╚';
        private const char BottomRight = '╝';
        private const char Horizontal = '═';
        private const char Vertical = '║';
        private const char VerticalRight = '╠';
        private const char VerticalLeft = '╣';
        private const char HorizontalDown = '╦';
        private const char HorizontalUp = '╩';
        private const char Cross = '╬';

        /// <summary>
        /// Creates a section header with consistent formatting
        /// </summary>
        public static string CreateHeader(string title, int width = 60)
        {
            if (string.IsNullOrEmpty(title))
                throw new ArgumentException("Title cannot be null or empty", nameof(title));

            title = title.Length > width - 4 ? title.Substring(0, width - 7) + "..." : title;
            int padding = (width - title.Length - 2) / 2;
            string paddedTitle = title.PadLeft(title.Length + padding).PadRight(width - 2);
            
            var sb = new StringBuilder();
            sb.AppendLine(new string(Horizontal, width));
            sb.AppendLine($"{Vertical}{paddedTitle}{Vertical}");
            sb.AppendLine(new string(Horizontal, width));
            return sb.ToString();
        }

        /// <summary>
        /// Creates a bordered box for content
        /// </summary>
        public static string CreateBox(string content, string? title = null, int minWidth = 60)
        {
            if (content == null) content = string.Empty;
            
            var lines = content.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            int width = Math.Max(minWidth, lines.Max(l => l.Length) + 4);
            
            var sb = new StringBuilder();
            
            // Top border
            if (!string.IsNullOrEmpty(title))
            {
                title = title.Length > width - 4 ? title.Substring(0, width - 7) + "..." : title;
                int padding = (width - title.Length - 2) / 2;
                string paddedTitle = title.PadLeft(title.Length + padding).PadRight(width - 2);
                
                sb.AppendLine($"{TopLeft}{new string(Horizontal, width - 2)}{TopRight}");
                sb.AppendLine($"{Vertical}{paddedTitle}{Vertical}");
                sb.AppendLine($"{VerticalRight}{new string(Horizontal, width - 2)}{VerticalLeft}");
            }
            else
            {
                sb.AppendLine($"{TopLeft}{new string(Horizontal, width - 2)}{TopRight}");
            }
            
            // Content
            foreach (var line in lines)
            {
                sb.AppendLine($"{Vertical} {line.PadRight(width - 4)} {Vertical}");
            }
            
            // Bottom border
            sb.AppendLine($"{BottomLeft}{new string(Horizontal, width - 2)}{BottomRight}");
            
            return sb.ToString();
        }

        /// <summary>
        /// Formats a key-value pair with consistent alignment
        /// </summary>
        public static string FormatKeyValue(string key, object value, int keyWidth = 25, int valueWidth = 15)
        {
            string keyStr = key?.ToString() ?? "";
            string valueStr = value?.ToString() ?? "";
            
            if (value is double doubleValue)
            {
                valueStr = doubleValue.ToString("F6").TrimEnd('0').TrimEnd('.');
                if (string.IsNullOrEmpty(valueStr)) valueStr = "0";
            }
            
            return $"{keyStr.PadRight(keyWidth)}: {valueStr.PadLeft(valueWidth)}";
        }

        /// <summary>
        /// Creates a table with consistent formatting
        /// </summary>
        public static string CreateTable(IEnumerable<string> headers, IEnumerable<IEnumerable<object>> rows)
        {
            var headerList = headers?.ToList() ?? new List<string>();
            var rowList = rows?.ToList() ?? new List<IEnumerable<object>>();
            
            // Convert all cells to strings and calculate column widths
            var stringRows = new List<List<string>>();
            var colWidths = new int[headerList.Count];
            
            // Process headers
            for (int i = 0; i < headerList.Count; i++)
            {
                colWidths[i] = Math.Max(colWidths[i], headerList[i]?.Length ?? 0);
            }
            
            // Process rows
            foreach (var row in rowList)
            {
                var stringRow = new List<string>();
                int col = 0;
                
                foreach (var cell in row)
                {
                    string cellStr = FormatCellValue(cell);
                    stringRow.Add(cellStr);
                    
                    if (col < colWidths.Length)
                    {
                        colWidths[col] = Math.Max(colWidths[col], cellStr.Length);
                    }
                    
                    col++;
                }
                
                stringRows.Add(stringRow);
            }
            
            // Build the table
            var sb = new StringBuilder();
            
            // Top border
            sb.Append(TopLeft);
            for (int i = 0; i < colWidths.Length; i++)
            {
                sb.Append(new string(Horizontal, colWidths[i] + 2));
                if (i < colWidths.Length - 1) sb.Append(HorizontalDown);
            }
            sb.AppendLine(TopRight.ToString());
            
            // Headers
            sb.Append(Vertical);
            for (int i = 0; i < headerList.Count; i++)
            {
                sb.Append(' ');
                sb.Append((headerList[i] ?? string.Empty).PadRight(colWidths[i]));
                sb.Append(' ');
                sb.Append(Vertical);
            }
            sb.AppendLine();
            
            // Header separator
            sb.Append(VerticalRight);
            for (int i = 0; i < colWidths.Length; i++)
            {
                sb.Append(new string(Horizontal, colWidths[i] + 2));
                if (i < colWidths.Length - 1) sb.Append(Cross);
            }
            sb.AppendLine(VerticalLeft.ToString());
            
            // Rows
            foreach (var row in stringRows)
            {
                sb.Append(Vertical);
                for (int i = 0; i < row.Count; i++)
                {
                    if (i < colWidths.Length)
                    {
                        sb.Append(' ');
                        sb.Append(row[i].PadRight(colWidths[i]));
                        sb.Append(' ');
                        sb.Append(Vertical);
                    }
                }
                sb.AppendLine();
            }
            
            // Bottom border
            sb.Append(BottomLeft);
            for (int i = 0; i < colWidths.Length; i++)
            {
                sb.Append(new string(Horizontal, colWidths[i] + 2));
                if (i < colWidths.Length - 1) sb.Append(HorizontalUp);
            }
            sb.AppendLine(BottomRight.ToString());
            
            return sb.ToString();
        }
        
        private static string FormatCellValue(object value)
        {
            if (value == null) return string.Empty;
            
            if (value is double d)
            {
                // Format doubles with 6 decimal places, but trim trailing zeros and decimal point if not needed
                return d.ToString("0.######");
            }
            
            return value.ToString();
        }
        
        /// <summary>
        /// Creates a horizontal rule with a consistent style
        /// </summary>
        public static string HorizontalRule(int width = 60, char lineChar = '─')
        {
            return new string(lineChar, width);
        }
        
        /// <summary>
        /// Formats a linear constraint into a readable string, e.g. "2 x1 + 3 x2 <= 5".
        /// Provides a single shared implementation so that all algorithms can call
        /// <c>OutputFormatter.FormatConstraint(...)</c>.
        /// </summary>
        public static string FormatConstraint(double[] coefficients, ConstraintType type, double rhs, string[]? variableNames = null, double epsilon = NumericalStabilityUtils.Epsilon)
        {
            if (coefficients == null) throw new ArgumentNullException(nameof(coefficients));

            var terms = new List<string>();
            for (int i = 0; i < coefficients.Length; i++)
            {
                double coef = coefficients[i];
                if (Math.Abs(coef) < epsilon) continue; // skip zero terms

                string varName = variableNames != null && i < variableNames.Length ? variableNames[i] : $"x{i + 1}";
                string sign = coef < 0 ? "- " : (terms.Count > 0 ? "+ " : string.Empty);
                double absCoef = Math.Abs(coef);
                string coefStr = absCoef.AlmostEquals(1.0, epsilon) ? string.Empty : absCoef.ToString("0.######");
                terms.Add($"{sign}{coefStr}{(coefStr != string.Empty ? " " : string.Empty)}{varName}".TrimStart());
            }
            if (terms.Count == 0) terms.Add("0");

            string rel = type switch
            {
                ConstraintType.LessThanOrEqual => "<=",
                ConstraintType.Equal => "=",
                ConstraintType.GreaterThanOrEqual => ">=",
                _ => "?"
            };

            return $"{string.Join(" ", terms)} {rel} {rhs.ToString("0.######")}";
        }

        /// <summary>
        /// Creates a section with a title and content
        /// </summary>
        public static string CreateSection(string title, string content)
        {
            var sb = new StringBuilder();
            sb.AppendLine(CreateHeader(title));
            sb.AppendLine(content);
            return sb.ToString();
        }
    }
}
