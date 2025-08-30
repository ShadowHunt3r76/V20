using System;
using System.Collections.Generic;
using System.Linq;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Represents a single iteration of the simplex algorithm with tableau information
    /// </summary>
    public class TableauIteration : ISimplexIteration
    {
        #region ISimplexIteration Implementation
        
        /// <inheritdoc/>
        public int Iteration { get; set; }
        
        /// <inheritdoc/>
        public double[,] Tableau { get; set; }
        
        /// <inheritdoc/>
        public int EnteringVariable { get; set; }
        
        /// <inheritdoc/>
        public int LeavingVariable { get; set; }
        
        /// <inheritdoc/>
        public double ObjectiveValue { get; set; }
        
        /// <inheritdoc/>
        public double[] Solution { get; set; }
        
        /// <inheritdoc/>
        public int[] Basis { get; set; }
        
        /// <inheritdoc/>
        public IList<string> BasisVariables { get; set; }
        
        /// <inheritdoc/>
        public string[] NonBasisVariables { get; set; }
        
        /// <inheritdoc/>
        public string[] VariableNames { get; set; }
        
        /// <inheritdoc/>
        public double[,] BasisInverse { get; set; }
        
        private bool _isOptimal;
        private bool _isUnbounded;
        private bool _isInfeasible;
        
        /// <inheritdoc/>
        public bool IsOptimal
        {
            get => _isOptimal || Status == "Optimal";
            set => _isOptimal = value;
        }
        
        /// <inheritdoc/>
        public bool IsUnbounded
        {
            get => _isUnbounded || Status == "Unbounded";
            set => _isUnbounded = value;
        }
        
        /// <inheritdoc/>
        public bool IsInfeasible
        {
            get => _isInfeasible || Status == "Infeasible";
            set => _isInfeasible = value;
        }
        
        /// <inheritdoc/>
        public double[] ReducedCosts { get; set; }
        
        /// <inheritdoc/>
        public double[] SimplexMultipliers { get; set; }
        
        /// <inheritdoc/>
        public string Description { get; set; }
        
        /// <inheritdoc/>
        public double PivotElement { get; set; }
        
        #endregion
        
        /// <summary>
        /// Gets or sets the status of the iteration
        /// </summary>
        public string Status { get; set; }
        
        /// <summary>
        /// Gets or sets the basic variable values
        /// </summary>
        public double[] BasicVariableValues { get; set; }
        
        /// <summary>
        /// Gets or sets the non-basic variable indices
        /// </summary>
        public int[] NonBasicVariableIndices { get; set; }
        
        private readonly List<int> _basisIndices = new List<int>();
        
        /// <summary>
        /// Gets or sets the list of basis indices
        /// </summary>
        public List<int> BasisIndices
        {
            get => _basisIndices;
            set
            {
                _basisIndices.Clear();
                if (value != null)
                {
                    _basisIndices.AddRange(value);
                }
            }
        }
        
        /// <summary>
        /// Initializes a new instance of the TableauIteration class
        /// </summary>
        public TableauIteration()
        {
            // Default constructor
        }
        
        /// <summary>
        /// Initializes a new instance of the TableauIteration class with the specified values
        /// </summary>
        public TableauIteration(
            int iteration,
            double[,] tableau,
            int enteringVariable,
            int leavingVariable,
            double objectiveValue,
            double[] solution,
            int[] basis,
            bool isOptimal = false,
            bool isUnbounded = false,
            bool isInfeasible = false)
        {
            Iteration = iteration;
            Tableau = (double[,])tableau.Clone();
            EnteringVariable = enteringVariable;
            LeavingVariable = leavingVariable;
            ObjectiveValue = objectiveValue;
            Solution = (double[])solution.Clone();
            Basis = (int[])basis.Clone();
            IsOptimal = isOptimal;
            IsUnbounded = isUnbounded;
            IsInfeasible = isInfeasible;
        }
    }
}
