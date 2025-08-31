using System;
using System.Collections.Generic;
using System.Linq;

namespace LinearProgramming.Parsing
{
    public class CanonicalLinearProgrammingModel
    {
        private double[][] _coefficientMatrix;
        private double[] _rhsVector;
        private double[] _objectiveCoefficients;
        private VariableType[] _variableTypes;
        private ConstraintType[] _constraintTypes;
        private string[] _variableNames;

        public double[][] CoefficientMatrix
        {
            get => _coefficientMatrix;
            set => _coefficientMatrix = value ?? throw new ArgumentNullException(nameof(value));
        }

        public double[] RHSVector
        {
            get => _rhsVector;
            set => _rhsVector = value ?? throw new ArgumentNullException(nameof(value));
        }

        public double[] RightHandSide
        {
            get => RHSVector;
            set => RHSVector = value;
        }

        public double[] ObjectiveCoefficients
        {
            get => _objectiveCoefficients;
            set => _objectiveCoefficients = value ?? throw new ArgumentNullException(nameof(value));
        }

        public VariableType[] VariableTypes
        {
            get => _variableTypes;
            set => _variableTypes = value ?? throw new ArgumentNullException(nameof(value));
        }

        public ConstraintType[] ConstraintTypes
        {
            get => _constraintTypes;
            set => _constraintTypes = value ?? throw new ArgumentNullException(nameof(value));
        }

        public string[] VariableNames
        {
            get => _variableNames;
            set => _variableNames = value ?? throw new ArgumentNullException(nameof(value));
        }

        public bool[] IsIntegerVariable => VariableTypes?.Select(t => t == VariableType.Integer || t == VariableType.Binary).ToArray() 
            ?? Array.Empty<bool>();

        public OptimizationType OptimizationType { get; set; } = OptimizationType.Maximize;
        public Dictionary<string, object> Metadata { get; set; } = new();

        public CanonicalLinearProgrammingModel()
        {
            _coefficientMatrix = Array.Empty<double[]>();
            _rhsVector = Array.Empty<double>();
            _objectiveCoefficients = Array.Empty<double>();
            _variableTypes = Array.Empty<VariableType>();
            _constraintTypes = Array.Empty<ConstraintType>();
            _variableNames = Array.Empty<string>();
        }
    }

    public enum OptimizationType
    {
        Minimize,
        Maximize
    }
}
