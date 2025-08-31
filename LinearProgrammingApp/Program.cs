using System;
using System.IO;
using System.Linq;
using System.Collections.Generic;
using LinearProgramming.Parsing;

namespace LinearProgrammingApp
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.Title = "Linear Programming Solver";
            Console.OutputEncoding = System.Text.Encoding.UTF8;
            
            try
            {
                // Get input and output file paths
                string inputFile = GetInputFile(args);
                if (inputFile == null) return;

                // Default output file is input filename with _output.txt suffix
                string outputFile = Path.Combine(
                    Path.GetDirectoryName(inputFile) ?? string.Empty,
                    $"{Path.GetFileNameWithoutExtension(inputFile)}_output.txt");

                // Parse the input file
                var parser = new ParsedLinearProgrammingModel.UniversalLinearProgrammingParser();
                var parsedModel = parser.ParseFromFile(inputFile);
                
                // Display welcome message
                Console.Clear();
                Console.WriteLine("=== Linear Programming Solver ===");
                Console.WriteLine($"Input file: {inputFile}");
                Console.WriteLine($"Output will be saved to: {outputFile}");
                Console.WriteLine($"Objective: {(parsedModel.Objective.Optimization == OptimizationType.Maximize ? "Maximize" : "Minimize")} {string.Join(" + ", parsedModel.Objective.Coefficients.Select((c, i) => $"{c}x{i+1}"))}");
                Console.WriteLine($"Variables: {parsedModel.Variables.Count} ({GetVariableTypeSummary(parsedModel)})");
                Console.WriteLine($"Constraints: {parsedModel.Constraints.Count}");
                
                // Start the main menu with the model and output file path
                SolverMenu.ShowMainMenu(parsedModel, outputFile);
            }
            catch (ParsedLinearProgrammingModel.LinearProgrammingParseException ex)
            {
                Console.WriteLine($"\nError parsing input file: {ex.Message}");
                if (ex.InnerException != null)
                {
                    Console.WriteLine($"Details: {ex.InnerException.Message}");
                }
                WaitForUser();
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nAn unexpected error occurred: {ex.Message}");
                Console.WriteLine($"\nStack trace: {ex.StackTrace}");
                WaitForUser();
            }
        }

        private static string GetInputFile(string[] args)
        {
            string inputFile = null;
            
            if (args.Length > 0)
            {
                // Check if first argument is an existing file
                if (File.Exists(args[0]))
                {
                    inputFile = Path.GetFullPath(args[0]);
                }
                // Check if first argument is a help flag
                else if (args[0] == "--help" || args[0] == "-h" || args[0] == "/?")
                {
                    Console.WriteLine("Linear Programming Solver");
                    Console.WriteLine("Usage: solve.exe [input_file] [output_file]");
                    Console.WriteLine("  input_file   Path to the input file containing the LP model");
                    Console.WriteLine("  output_file  (Optional) Path to save the output (default: <input_filename>_output.txt)");
                    Console.WriteLine("\nIf no arguments are provided, you will be prompted for the input file.");
                    return null;
                }
            }

            // If no valid file was provided as argument, prompt user
            while (string.IsNullOrEmpty(inputFile) || !File.Exists(inputFile))
            {
                if (!string.IsNullOrEmpty(inputFile))
                {
                    Console.WriteLine($"File not found: {inputFile}");
                }
                
                Console.Write("\nEnter the path to the input file (or press Enter to exit): ");
                inputFile = Console.ReadLine()?.Trim('"');
                
                if (string.IsNullOrEmpty(inputFile))
                {
                    Console.WriteLine("No input file specified. Exiting...");
                    return null;
                }
                
                if (File.Exists(inputFile))
                {
                    break;
                }
            }
            
            return Path.GetFullPath(inputFile).Trim('"');
        }

        private static string GetOutputFile(string[] args, string inputFile)
        {
            if (args.Length > 1)
            {
                return Path.GetFullPath(args[1]);
            }
            else
            {
                return Path.Combine(
                    Path.GetDirectoryName(inputFile) ?? string.Empty,
                    $"{Path.GetFileNameWithoutExtension(inputFile)}_output.txt");
            }
        }
        
        private static string GetVariableTypeSummary(ParsedLinearProgrammingModel model)
        {
            // Since the current VariableType enum doesn't distinguish between integer and binary variables,
            // we'll just count the total number of variables
            int totalVariables = model.Variables.Count;
            
            return $"{totalVariables} variables";
        }
        
        private static void WaitForUser()
        {
            Console.WriteLine("\nPress any key to exit...");
            Console.ReadKey();
        }
    }
}
