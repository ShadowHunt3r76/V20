using System;
using System.Text;

namespace LinearProgramming.Algorithms.SensitivityAnalysis
{
    /// <summary>
    /// Provides common output formatting for sensitivity analysis results
    /// </summary>
    public static class SensitivityOutputFormatter
    {
        public static void PrintSectionHeader(string title, char separator = '=', int width = 80)
        {
            if (string.IsNullOrEmpty(title))
                throw new ArgumentNullException(nameof(title));
                
            string separatorLine = new string(separator, width);
            Console.WriteLine("\n" + separatorLine);
            Console.WriteLine(title);
            Console.WriteLine(separatorLine);
        }

        public static void PrintSubsectionHeader(string title, char separator = '-', int width = 60)
        {
            if (string.IsNullOrEmpty(title))
                throw new ArgumentNullException(nameof(title));
                
            string separatorLine = new string(separator, Math.Min(title.Length + 4, width));
            Console.WriteLine("\n" + separatorLine);
            Console.WriteLine($"  {title}  ");
            Console.WriteLine(separatorLine);
        }

        public static string FormatTableRow(params (string, int)[] columns)
        {
            var sb = new StringBuilder();
            foreach (var (text, width) in columns)
            {
                string padded = (text ?? string.Empty).PadRight(width);
                if (padded.Length > width)
                    padded = padded.Substring(0, width - 3) + "...";
                sb.Append(padded);
            }
            return sb.ToString();
        }

        public static string FormatKeyValue(string key, object value, int keyWidth = 30, int valueWidth = 30)
        {
            return $"{key.PadRight(keyWidth)}: {value}";
        }

        public static string FormatRange(double? lower, double? upper, string format = "F6")
        {
            if (!lower.HasValue && !upper.HasValue)
                return "Unbounded";
                
            string lowerStr = lower.HasValue 
                ? lower.Value.ToString(format) 
                : "-∞";
                
            string upperStr = upper.HasValue 
                ? upper.Value.ToString(format) 
                : "+∞";
                
            return $"[{lowerStr}, {upperStr}]";
        }
    }
}
