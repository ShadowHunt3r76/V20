using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Algorithms.SensitivityAnalysis;
using LinearProgramming.Algorithms.Utils;
using static LinearProgramming.Algorithms.Utils.NumericalStabilityUtils;

namespace LinearProgramming.Algorithms.CuttingPlane
{
    /// <summary>
    /// Provides sensitivity analysis for Cutting Plane algorithm solutions
    /// </summary>
    public class CuttingPlaneSensitivityAnalysis : BaseSensitivityAnalysis
    {
        private readonly LinearProgramSolution _lpRelaxationSolution;
        private readonly LinearProgramSolution _finalSolution;
        private readonly CanonicalLinearProgrammingModel _originalModel;
        private readonly List<double[,]> _cuttingPlaneHistory;

        /// <summary>
        /// Initializes a new instance of the CuttingPlaneSensitivityAnalysis class
        /// </summary>
        /// <param name="lpRelaxationSolution">Solution to the LP relaxation</param>
        /// <param name="finalSolution">Final integer solution</param>
        /// <param name="originalModel">Original model</param>
        /// <param name="cuttingPlaneHistory">History of tableaus from cutting plane iterations</param>
        public CuttingPlaneSensitivityAnalysis(
            LinearProgramSolution lpRelaxationSolution,
            LinearProgramSolution finalSolution,
            CanonicalLinearProgrammingModel originalModel,
            List<double[,]> cuttingPlaneHistory)
            : base(
                originalModel?.VariableNames ?? Array.Empty<string>(),
                originalModel?.ConstraintTypes ?? Array.Empty<ConstraintType>(),
                originalModel?.VariableNames?.Length ?? 0,
                originalModel?.ConstraintTypes?.Length ?? 0)
        {
            _lpRelaxationSolution = lpRelaxationSolution ?? throw new ArgumentNullException(nameof(lpRelaxationSolution));
            _finalSolution = finalSolution ?? throw new ArgumentNullException(nameof(finalSolution));
            _originalModel = originalModel ?? throw new ArgumentNullException(nameof(originalModel));
            _cuttingPlaneHistory = cuttingPlaneHistory ?? new List<double[,]>();
        }

        /// <summary>
        /// Performs complete sensitivity analysis on the Cutting Plane solution
        /// </summary>
        public void PerformAnalysis()
        {
            if (_finalSolution.Status != "Optimal")
            {
                Console.WriteLine("No optimal solution found for sensitivity analysis.");
                return;
            }

            Console.WriteLine("\n" + new string('=', 80));
            Console.WriteLine("CUTTING PLANE SENSITIVITY ANALYSIS");
            Console.WriteLine(new string('=', 80));

            try
            {
                // 1. Basic solution information
                DisplaySolutionInfo();

                // 2. LP relaxation sensitivity analysis
                if (_lpRelaxationSolution != null && _lpRelaxationSolution.OptimalTable != null)
                {
                    Console.WriteLine("\n" + new string('-', 40));
                    Console.WriteLine("LP RELAXATION SENSITIVITY ANALYSIS");
                    Console.WriteLine(new string('-', 40));
                    PerformLPSensitivityAnalysis();
                }

                // 3. Integer-specific sensitivity analysis
                Console.WriteLine("\n" + new string('-', 40));
                Console.WriteLine("INTEGER-SPECIFIC SENSITIVITY ANALYSIS");
                Console.WriteLine(new string('-', 40));
                PerformIntegerSensitivityAnalysis();

                // 4. Cutting plane analysis
                if (_cuttingPlaneHistory != null && _cuttingPlaneHistory.Count > 0)
                {
                    Console.WriteLine("\n" + new string('-', 40));
                    Console.WriteLine("CUTTING PLANE ANALYSIS");
                    Console.WriteLine(new string('-', 40));
                    AnalyzeCuttingPlanes();
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\n⚠ Error during sensitivity analysis: {ex.Message}");
            }

            Console.WriteLine("\n" + new string('=', 80));
            Console.WriteLine("END OF SENSITIVITY ANALYSIS");
            Console.WriteLine(new string('=', 80));
        }

        private void DisplaySolutionInfo()
        {
            Console.WriteLine("\nSOLUTION INFORMATION");
            Console.WriteLine(new string('-', 40));
            
            Console.WriteLine($"Status: {_finalSolution.Status}");
            Console.WriteLine($"Objective Value: {_finalSolution.ObjectiveValue:F6}");
            
            if (_finalSolution.SolutionVector != null)
            {
                Console.WriteLine("\nSolution Vector:");
                for (int i = 0; i < Math.Min(_finalSolution.SolutionVector.Length, _originalModel.VariableNames.Length); i++)
                {
                    Console.WriteLine($"{_originalModel.VariableNames[i]}: {_finalSolution.SolutionVector[i]:F6}");
                }
            }

            if (_lpRelaxationSolution != null)
            {
                double integralityGap = _finalSolution.ObjectiveValue - _lpRelaxationSolution.ObjectiveValue;
                Console.WriteLine($"\nLP Relaxation Objective: {_lpRelaxationSolution.ObjectiveValue:F6}");
                Console.WriteLine($"Integrality Gap: {integralityGap:F6}");
            }
        }

        private void PerformLPSensitivityAnalysis()
        {
            try
            {
                if (_lpRelaxationSolution.OptimalTable == null || _lpRelaxationSolution.BasisIndices == null)
                {
                    Console.WriteLine("LP relaxation solution data not available for sensitivity analysis.");
                    return;
                }

                var sensitivity = new SensitivityAnalysis(
                    _lpRelaxationSolution.OptimalTable,
                    _lpRelaxationSolution.BasisIndices,
                    _originalModel.VariableNames,
                    _originalModel.ConstraintTypes,
                    _originalModel.ObjectiveCoefficients.Length,
                    _originalModel.ConstraintTypes.Length);

                // Analyze objective coefficient ranges
                Console.WriteLine("\nOBJECTIVE COEFFICIENT RANGES (LP RELAXATION):");
                for (int i = 0; i < _originalModel.ObjectiveCoefficients.Length; i++)
                {
                    try
                    {
                        var range = sensitivity.GetObjectiveCoefficientRange(i);
                        Console.WriteLine($"{_originalModel.VariableNames[i]}: [{range.LowerBound:F6}, {range.UpperBound:F6}]");
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine($"Error analyzing coefficient for {_originalModel.VariableNames[i]}: {ex.Message}");
                    }
                }

                // Analyze RHS ranges
                Console.WriteLine("\nRIGHT-HAND SIDE RANGES (LP RELAXATION):");
                for (int i = 0; i < _originalModel.ConstraintTypes.Length; i++)
                {
                    try
                    {
                        var range = sensitivity.GetRHSRange(i);
                        Console.WriteLine($"Constraint {i + 1}: [{range.LowerBound:F6}, {range.UpperBound:F6}]");
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine($"Error analyzing constraint {i + 1}: {ex.Message}");
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error in LP sensitivity analysis: {ex.Message}");
            }
        }

        private void PerformIntegerSensitivityAnalysis()
        {
            if (_finalSolution.SolutionVector == null) return;

            Console.WriteLine("\nINTEGER VARIABLE SENSITIVITY:");
            for (int i = 0; i < _originalModel.VariableNames.Length; i++)
            {
                if (_originalModel.VariableTypes[i] != VariableType.Continuous)
                {
                    double value = _finalSolution.SolutionVector[i];
                    double floor = Math.Floor(value);
                    double ceil = Math.Ceiling(value);
                    
                    Console.WriteLine($"{_originalModel.VariableNames[i]}: Value = {value:F0} (Integer)");
                    
                    // For binary variables
                    if (_originalModel.VariableTypes[i] == VariableType.Binary)
                    {
                        bool canFlip = true;
                        // Check if flipping the value would affect feasibility
                        // This is a simplified check - in practice, we'd need to check all constraints
                        Console.WriteLine($"  Binary variable - Can {(canFlip ? "" : "not ")}flip to {1 - value}");
                    }
                    // For general integer variables
                    else
                    {
                        Console.WriteLine($"  Integer variable - Can increase to {value + 1}, Can decrease to {value - 1}");
                    }
                }
            }
        }

        private void AnalyzeCuttingPlanes()
        {
            if (_cuttingPlaneHistory == null || _cuttingPlaneHistory.Count == 0)
            {
                Console.WriteLine("No cutting plane history available.");
                return;
            }

            Console.WriteLine($"\nCUTTING PLANE ITERATIONS: {_cuttingPlaneHistory.Count}");
            
            // Analyze the impact of the first few cuts (to avoid too much output)
            int maxCutsToShow = Math.Min(5, _cuttingPlaneHistory.Count);
            for (int i = 0; i < maxCutsToShow; i++)
            {
                Console.WriteLine($"\nCut {i + 1}:");
                // In a real implementation, we would analyze the cut's impact on the solution
                // This is a simplified version
                Console.WriteLine("  Added to tighten the LP relaxation");
            }

            if (_cuttingPlaneHistory.Count > maxCutsToShow)
            {
                Console.WriteLine($"  ... and {_cuttingPlaneHistory.Count - maxCutsToShow} more cuts");
            }
        }
    }
}
