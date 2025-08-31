using System;
using System.Collections.Generic;
using System.Linq;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Represents an iteration of the revised simplex method
    /// </summary>
    public class RevisedSimplexIteration : ISimplexIteration
    {
        private string _status;
        private readonly double[,] _tableau;
        private bool _isOptimal;
        private bool _isUnbounded;
        private bool _isInfeasible;
        private readonly List<int> _basisIndices = new List<int>();
        
        /// <summary>
        /// Initializes a new instance of the <see cref="RevisedSimplexIteration"/> class
        /// </summary>
        public RevisedSimplexIteration()
        {
            _tableau = new double[0, 0]; // Empty tableau for revised simplex
            BasisVariables = new List<string>();
            NonBasisVariables = Array.Empty<string>();
            VariableNames = Array.Empty<string>();
            Description = "";
        }
        
        /// <summary>
        /// Gets the iteration number
        /// </summary>
        public int Iteration { get; set; }
        
        /// <summary>
        /// Gets the current tableau (returns empty array for revised simplex)
        /// </summary>
        public double[,] Tableau => _tableau;
        
        /// <summary>
        /// Gets the index of the entering variable (-1 if none)
        /// </summary>
        public int EnteringVariable { get; set; } = -1;
        
        /// <summary>
        /// Gets the index of the leaving variable (-1 if none)
        /// </summary>
        public int LeavingVariable { get; set; } = -1;
        
        /// <summary>
        /// Gets the current objective value
        /// </summary>
        public double ObjectiveValue { get; set; }
        
        /// <summary>
        /// Gets the current solution vector
        /// </summary>
        public double[] Solution { get; set; } = Array.Empty<double>();
        
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
        /// Gets the current basis variables as an array
        /// </summary>
        public int[] Basis => _basisIndices.ToArray();
        
        /// <summary>
        /// Gets the names of the basis variables
        /// </summary>
        public IList<string> BasisVariables { get; set; }
        
        /// <summary>
        /// Gets the names of the non-basis variables
        /// </summary>
        public string[] NonBasisVariables { get; set; }
        
        /// <summary>
        /// Gets the names of all variables
        /// </summary>
        public string[] VariableNames { get; set; }
        
        /// <summary>
        /// Gets the description of the iteration
        /// </summary>
        public string Description { get; set; }
        
        /// <summary>
        /// Gets the pivot element used in this iteration
        /// </summary>
        public double PivotElement { get; set; }
        
        /// <summary>
        /// Gets or sets the status of the iteration
        /// </summary>
        public string Status 
        { 
            get => _status;
            set 
            {
                _status = value ?? throw new ArgumentNullException(nameof(value));
                // Update status flags based on status string
                _isOptimal = _status == "Optimal";
                _isUnbounded = _status == "Unbounded";
                _isInfeasible = _status == "Infeasible";
            }
        }
        
        /// <summary>
        /// Gets whether the solution is optimal
        /// </summary>
        public bool IsOptimal => _isOptimal;
        
        /// <summary>
        /// Gets whether the problem is unbounded
        /// </summary>
        public bool IsUnbounded => _isUnbounded;
        
        /// <summary>
        /// Gets whether the problem is infeasible
        /// </summary>
        public bool IsInfeasible => _isInfeasible;
        
        /// <summary>
        /// Gets or sets the reduced costs
        /// </summary>
        public double[] ReducedCosts { get; set; } = Array.Empty<double>();
        
        /// <summary>
        /// Gets or sets the simplex multipliers (dual variables)
        /// </summary>
        public double[] SimplexMultipliers { get; set; } = Array.Empty<double>();
        
        /// <summary>
        /// Gets or sets the inverse of the basis matrix
        /// </summary>
        public double[,] BasisInverse { get; set; } = new double[0, 0];
        
        /// <summary>
        /// Gets or sets the pivot value used in this iteration
        /// </summary>
        public double PivotValue { get; set; }
        
        // Additional properties for revised simplex method
        public string Type { get; set; }
        public double[] BasicSolution { get; set; } = Array.Empty<double>();
    }
}
