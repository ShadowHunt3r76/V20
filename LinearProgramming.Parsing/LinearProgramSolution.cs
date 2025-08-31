using System;
using System.Collections.Generic;

namespace LinearProgramming.Parsing
{
    public class LinearProgramSolution : ILinearProgramSolution
    {
        public string Status { get; set; } = string.Empty;
        public double ObjectiveValue { get; set; }
        public Dictionary<string, double> VariableValues { get; set; } = new();
        public double[] SolutionVector { get; set; } = Array.Empty<double>();
        public double[] DualValues { get; set; } = Array.Empty<double>();
        public double[] ReducedCosts { get; set; } = Array.Empty<double>();
        public int Iterations { get; set; }
        public TimeSpan SolveTime { get; set; }
        public string Algorithm { get; set; } = string.Empty;
        public string Message { get; set; } = string.Empty;
        
        // ILinearProgramSolution interface implementation
        public double[,] InitialTable { get; set; } = new double[0, 0];
        public double[,] OptimalTable { get; set; } = new double[0, 0];
        public List<double[,]> TableHistory { get; set; } = new();
        public string[] BasisVariables { get; set; } = Array.Empty<string>();
        public int[] BasisIndices { get; set; } = Array.Empty<int>();
        public double[,] CanonicalMatrix { get; set; } = new double[0, 0];
        public string[] VariableNames { get; set; } = Array.Empty<string>();
    }
}
