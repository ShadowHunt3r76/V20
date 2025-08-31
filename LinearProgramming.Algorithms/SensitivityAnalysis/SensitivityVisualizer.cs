using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LinearProgramming.Algorithms.SensitivityAnalysis
{
    /// <summary>
    /// Provides visualization capabilities for sensitivity analysis results
    /// </summary>
    public static class SensitivityVisualizer
    {
        /// <summary>
        /// Creates a text-based bar chart for visualizing value ranges
        /// </summary>
        public static string CreateBarChart(string title, Dictionary<string, (double min, double current, double max)> data, int width = 50)
        {
            if (data == null || data.Count == 0)
                return "No data available for visualization";

            var sb = new StringBuilder();
            sb.AppendLine($"{title}\n{new string('=', title.Length)}");

            // Find the maximum length of the variable names for alignment
            int maxNameLength = data.Keys.Max(n => n.Length);
            
            foreach (var item in data)
            {
                string name = item.Key.PadRight(maxNameLength);
                double min = item.Value.min;
                double current = item.Value.current;
                double max = item.Value.max;

                // Skip if the range is not meaningful
                if (double.IsNegativeInfinity(min) && double.IsPositiveInfinity(max))
                {
                    sb.AppendLine($"{name}: No effective bounds");
                    continue;
                }

                // Calculate positions for the current value and range markers
                double range = max - min;
                int currentPos = width / 2; // Default position if range is 0
                
                if (range > 0)
                {
                    currentPos = (int)((current - min) / range * (width - 1));
                    currentPos = Math.Max(0, Math.Min(width - 1, currentPos));
                }

                // Create the bar
                char[] bar = new char[width];
                for (int i = 0; i < width; i++)
                {
                    if (i == currentPos)
                        bar[i] = '|';  // Current value marker
                    else if (i == 0)
                        bar[i] = '[';   // Start of range
                    else if (i == width - 1)
                        bar[i] = ']';   // End of range
                    else if (i < currentPos)
                        bar[i] = '=';   // Before current value
                    else
                        bar[i] = '-';   // After current value
                }

                // Add the formatted line
                sb.AppendLine($"{name}: [{min,8:F2}, {max,8:F2}] {new string(bar)} {current,8:F2}");
            }

            return sb.ToString();
        }

        /// <summary>
        /// Creates a text-based histogram for visualizing value distributions
        /// </summary>
        public static string CreateHistogram(string title, Dictionary<string, double> data, int height = 10, int width = 50)
        {
            if (data == null || data.Count == 0)
                return "No data available for histogram";

            var sb = new StringBuilder();
            sb.AppendLine($"{title}\n{new string('=', title.Length)}");

            // Find the maximum value for scaling
            double maxValue = data.Values.Max();
            if (maxValue <= 0) maxValue = 1; // Avoid division by zero

            // Sort data by value (descending)
            var sortedData = data.OrderByDescending(kv => kv.Value).ToList();
            
            // Find the maximum key length for alignment
            int maxKeyLength = sortedData.Max(kv => kv.Key.Length);

            foreach (var item in sortedData)
            {
                string key = item.Key.PadRight(maxKeyLength);
                double value = item.Value;
                int barLength = (int)((value / maxValue) * width);
                
                sb.AppendLine($"{key} |{new string('█', barLength)} {value,8:F4}");
            }

            return sb.ToString();
        }

        /// <summary>
        /// Creates a comparison table for before/after analysis
        /// </summary>
        public static string CreateComparisonTable(
            string title, 
            Dictionary<string, (double original, double newValue, string unit)> data,
            string[]? columns = null)
        {
            if (data == null || data.Count == 0)
                return "No data available for comparison";

            var columnNames = columns ?? new[] { "Variable", "Original", "New", "Change", "% Change" };
            
            var sb = new StringBuilder();
            sb.AppendLine($"{title}\n{new string('=', title.Length)}");

            // Calculate column widths
            int[] colWidths = new int[columnNames.Length];
            for (int i = 0; i < columnNames.Length; i++)
                colWidths[i] = columnNames[i].Length;

            // Update column widths based on data
            foreach (var item in data)
            {
                colWidths[0] = Math.Max(colWidths[0], item.Key.Length);
                colWidths[1] = Math.Max(colWidths[1], item.Value.original.ToString("F4").Length);
                colWidths[2] = Math.Max(colWidths[2], item.Value.newValue.ToString("F4").Length);
                
                double change = item.Value.newValue - item.Value.original;
                double pctChange = item.Value.original != 0 ? (change / Math.Abs(item.Value.original)) * 100 : 0;
                
                colWidths[3] = Math.Max(colWidths[3], change.ToString("+0.00;-0.00;0").Length);
                colWidths[4] = Math.Max(colWidths[4], pctChange.ToString("+0.00%;-0.00%;0%").Length);
            }

            // Add some padding
            for (int i = 0; i < colWidths.Length; i++)
                colWidths[i] += 2;

            // Print header
            for (int i = 0; i < columns.Length; i++)
                sb.Append(columns[i].PadRight(colWidths[i]));
            sb.AppendLine();
            sb.AppendLine(new string('-', colWidths.Sum()));

            // Print rows
            foreach (var item in data)
            {
                double change = item.Value.newValue - item.Value.original;
                double pctChange = item.Value.original != 0 ? (change / Math.Abs(item.Value.original)) * 100 : 0;
                
                sb.Append(item.Key.PadRight(colWidths[0]));
                sb.Append(item.Value.original.ToString("F4").PadRight(colWidths[1]));
                sb.Append(item.Value.newValue.ToString("F4").PadRight(colWidths[2]));
                sb.Append(change.ToString("+0.00;-0.00;0").PadRight(colWidths[3]));
                sb.AppendLine(pctChange.ToString("+0.00%;-0.00%;0%").PadRight(colWidths[4]));
            }

            return sb.ToString();
        }
    }
}
