using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.PrimalSimplex;
using static LinearProgramming.Algorithms.PrimalSimplex.MatrixUtils;

namespace LinearProgramming.Algorithms.BranchAndBound
{
    using CanonicalLinearProgrammingModel = LinearProgramming.Parsing.ParsedLinearProgrammingModel.CanonicalLinearProgrammingModel;
    using LinearProgramSolution = LinearProgramming.Algorithms.PrimalSimplex.LinearProgramSolution;

    /// <summary>
    /// Represents a node in the branch and bound tree
    /// </summary>
    public class BranchAndBoundNode
    {
        public int Id { get; set; }
        public int? ParentId { get; set; }
        public CanonicalLinearProgrammingModel Model { get; set; }
        public LinearProgramSolution Solution { get; set; }
        public string BranchingConstraint { get; set; }
        public double UpperBound { get; set; }
        public bool IsFathomed { get; set; }
        public string FathomReason { get; set; }
        public int BranchedVariableIndex { get; set; }
        public List<double[,]> TableIterations { get; set; } = new List<double[,]>();
        public int Depth { get; set; }
    }

    /// <summary>
    /// Result class containing all branch and bound solving information
    /// </summary>
    public class BranchAndBoundSolution
    {
        public List<BranchAndBoundNode> AllNodes { get; set; } = new List<BranchAndBoundNode>();
        public BranchAndBoundNode BestCandidate { get; set; }
        public double BestObjectiveValue { get; set; } = double.MinValue;
        public double[] BestSolution { get; set; }
        public string OverallStatus { get; set; }
        public CanonicalLinearProgrammingModel CanonicalForm { get; set; }
        public double[,] CanonicalMatrix { get; set; }
        public int TotalNodesExplored { get; set; }
        public int TotalNodesFathomed { get; set; }
    }

    /// <summary>
    /// Branch and Bound Simplex Algorithm solver for Integer Programming problems
    /// Implements backtracking, complete sub-problem generation, node fathoming, and table iteration tracking
    /// </summary>
    public class BranchAndBoundSolver
    {
        private const double EPSILON = 1e-6;
        private int nextNodeId = 1;
        private PrimalSimplexSolver simplexSolver = new PrimalSimplexSolver();

        /// <summary>
        /// Solves the integer programming problem using Branch and Bound with Simplex Algorithm
        /// </summary>
        /// <param name="originalModel">The parsed linear programming model with integer/binary variables</param>
        /// <returns>Complete branch and bound solution with all sub-problems and iterations</returns>
        public BranchAndBoundSolution Solve(ParsedLinearProgrammingModel originalModel)
        {
            var solution = new BranchAndBoundSolution();
            
            // Convert to canonical form and store for display
            solution.CanonicalForm = originalModel.ToCanonicalForm();
            solution.CanonicalMatrix = MatrixUtils.ConvertToMatrix(solution.CanonicalForm.CoefficientMatrix);
            
            // Check if any variables are integer or binary
            bool hasIntegerVariables = originalModel.Variables.Any(v => v.Type == VariableType.Integer || v.Type == VariableType.Binary);
            if (!hasIntegerVariables)
            {
                // If no integer variables, solve as regular LP
                var lpSolution = simplexSolver.Solve(solution.CanonicalForm);
                solution.BestCandidate = new BranchAndBoundNode
                {
                    Id = 0,
                    Model = solution.CanonicalForm,
                    Solution = lpSolution
                };
                solution.BestObjectiveValue = lpSolution.ObjectiveValue;
                solution.BestSolution = lpSolution.SolutionVector;
                solution.OverallStatus = lpSolution.Status;
                return solution;
            }

            // Initialize with root node
            var rootNode = CreateRootNode(solution.CanonicalForm);
            solution.AllNodes.Add(rootNode);
            
            // Stack for backtracking (DFS approach)
            var nodeStack = new Stack<BranchAndBoundNode>();
            nodeStack.Push(rootNode);

            solution.TotalNodesExplored = 0;
            solution.TotalNodesFathomed = 0;

            // Branch and Bound main loop with backtracking
            while (nodeStack.Count > 0)
            {
                var currentNode = nodeStack.Pop();
                solution.TotalNodesExplored++;

                Console.WriteLine($"\n=== Processing Node {currentNode.Id} (Depth {currentNode.Depth}) ===");
                if (currentNode.BranchingConstraint != null)
                {
                    Console.WriteLine($"Branching Constraint: {currentNode.BranchingConstraint}");
                }

                // Solve the current subproblem using simplex
                currentNode.Solution = simplexSolver.Solve(currentNode.Model);
                
                // Store all table iterations for this node
                currentNode.TableIterations.AddRange(currentNode.Solution.TableHistory);
                if (currentNode.Solution.InitialTable != null)
                {
                    currentNode.TableIterations.Insert(0, currentNode.Solution.InitialTable);
                }
                if (currentNode.Solution.OptimalTable != null)
                {
                    currentNode.TableIterations.Add(currentNode.Solution.OptimalTable);
                }

                Console.WriteLine($"Node {currentNode.Id} Status: {currentNode.Solution.Status}");
                if (currentNode.Solution.Status == "Optimal")
                {
                    Console.WriteLine($"Objective Value: {currentNode.Solution.ObjectiveValue:F3}");
                }

                // Check fathoming conditions
                if (ShouldFathomNode(currentNode, solution.BestObjectiveValue, originalModel))
                {
                    currentNode.IsFathomed = true;
                    solution.TotalNodesFathomed++;
                    Console.WriteLine($"Node {currentNode.Id} FATHOMED: {currentNode.FathomReason}");
                    continue;
                }

                // Check if solution is integer-feasible
                if (IsIntegerFeasible(currentNode.Solution.SolutionVector, originalModel))
                {
                    // Update best solution if this is better
                    if (currentNode.Solution.ObjectiveValue > solution.BestObjectiveValue)
                    {
                        solution.BestObjectiveValue = currentNode.Solution.ObjectiveValue;
                        solution.BestSolution = (double[])currentNode.Solution.SolutionVector.Clone();
                        solution.BestCandidate = currentNode;
                        Console.WriteLine($"*** NEW BEST INTEGER SOLUTION FOUND ***");
                        Console.WriteLine($"Objective Value: {solution.BestObjectiveValue:F3}");
                        PrintSolution(solution.BestSolution);
                    }
                    // Fathom by integrality
                    currentNode.IsFathomed = true;
                    currentNode.FathomReason = "Integer solution found";
                    solution.TotalNodesFathomed++;
                    continue;
                }

                // Branch on fractional integer/binary variable
                var branchingVariable = FindBranchingVariable(currentNode.Solution.SolutionVector, originalModel);
                if (branchingVariable == -1)
                {
                    // This should not happen if IsIntegerFeasible returned false
                    currentNode.IsFathomed = true;
                    currentNode.FathomReason = "No fractional integer variables found";
                    solution.TotalNodesFathomed++;
                    continue;
                }

                // Create child nodes (backtracking will handle them via stack)
                var childNodes = CreateChildNodes(currentNode, branchingVariable, originalModel);
                
                // Add children to stack (reverse order for consistent left-first exploration)
                for (int i = childNodes.Count - 1; i >= 0; i--)
                {
                    solution.AllNodes.Add(childNodes[i]);
                    nodeStack.Push(childNodes[i]);
                }

                Console.WriteLine($"Created {childNodes.Count} child nodes from Node {currentNode.Id}");
            }

            // Determine overall status
            if (solution.BestCandidate != null)
            {
                solution.OverallStatus = "Optimal Integer Solution Found";
            }
            else if (solution.AllNodes.Any(n => n.Solution?.Status == "Unbounded"))
            {
                solution.OverallStatus = "Integer Problem is Unbounded";
            }
            else
            {
                solution.OverallStatus = "No Integer Solution Exists";
            }

            Console.WriteLine($"\n=== BRANCH AND BOUND COMPLETE ===");
            Console.WriteLine($"Total Nodes Explored: {solution.TotalNodesExplored}");
            Console.WriteLine($"Total Nodes Fathomed: {solution.TotalNodesFathomed}");
            Console.WriteLine($"Overall Status: {solution.OverallStatus}");

            return solution;
        }

        /// <summary>
        /// Creates the root node of the branch and bound tree
        /// </summary>
        private BranchAndBoundNode CreateRootNode(CanonicalLinearProgrammingModel canonicalModel)
        {
            return new BranchAndBoundNode
            {
                Id = 0,
                ParentId = null,
                Model = CloneModel(canonicalModel),
                BranchingConstraint = "Root Node (Original Problem)",
                Depth = 0
            };
        }

        /// <summary>
        /// Determines if a node should be fathomed based on bound, infeasibility, or unboundedness
        /// </summary>
        private bool ShouldFathomNode(BranchAndBoundNode node, double bestKnownValue, ParsedLinearProgrammingModel originalModel)
        {
            // Fathom by infeasibility
            if (node.Solution.Status == "Infeasible")
            {
                node.FathomReason = "Infeasible";
                return true;
            }

            // Fathom by unboundedness (for maximization problems, this means no optimal integer solution in this branch)
            if (node.Solution.Status == "Unbounded")
            {
                node.FathomReason = "Unbounded (LP relaxation)";
                return true;
            }

            // Fathom by bound (for maximization: if upper bound <= best known integer value)
            if (node.Solution.Status == "Optimal")
            {
                node.UpperBound = node.Solution.ObjectiveValue;
                if (bestKnownValue != double.MinValue && node.UpperBound <= bestKnownValue + EPSILON)
                {
                    node.FathomReason = $"Bound ({node.UpperBound:F3} ? {bestKnownValue:F3})";
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Checks if the solution satisfies integer constraints
        /// </summary>
        private bool IsIntegerFeasible(double[] solution, ParsedLinearProgrammingModel originalModel)
        {
            if (solution == null) return false;

            for (int i = 0; i < Math.Min(solution.Length, originalModel.Variables.Count); i++)
            {
                var varType = originalModel.Variables[i].Type;
                
                if (varType == VariableType.Integer)
                {
                    if (Math.Abs(solution[i] - Math.Round(solution[i])) > EPSILON)
                        return false;
                }
                else if (varType == VariableType.Binary)
                {
                    double rounded = Math.Round(solution[i]);
                    if (Math.Abs(solution[i] - rounded) > EPSILON || (rounded != 0 && rounded != 1))
                        return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Finds the first fractional integer or binary variable to branch on
        /// </summary>
        private int FindBranchingVariable(double[] solution, ParsedLinearProgrammingModel originalModel)
        {
            for (int i = 0; i < Math.Min(solution.Length, originalModel.Variables.Count); i++)
            {
                var varType = originalModel.Variables[i].Type;
                
                if (varType == VariableType.Integer || varType == VariableType.Binary)
                {
                    if (Math.Abs(solution[i] - Math.Round(solution[i])) > EPSILON)
                        return i;
                }
            }
            return -1;
        }

        /// <summary>
        /// Creates child nodes by adding branching constraints
        /// </summary>
        private List<BranchAndBoundNode> CreateChildNodes(BranchAndBoundNode parentNode, int variableIndex, ParsedLinearProgrammingModel originalModel)
        {
            var children = new List<BranchAndBoundNode>();
            double fractionalValue = parentNode.Solution.SolutionVector[variableIndex];
            var varType = originalModel.Variables[variableIndex].Type;

            if (varType == VariableType.Binary)
            {
                // For binary variables, create x ? 0 and x ? 1 branches
                children.Add(CreateChildWithConstraint(parentNode, variableIndex, 0, true, "? 0"));  // x ? 0
                children.Add(CreateChildWithConstraint(parentNode, variableIndex, 1, false, "? 1")); // x ? 1
            }
            else if (varType == VariableType.Integer)
            {
                // For integer variables, create x ? floor(value) and x ? ceil(value) branches
                int floorValue = (int)Math.Floor(fractionalValue);
                int ceilValue = (int)Math.Ceiling(fractionalValue);
                
                children.Add(CreateChildWithConstraint(parentNode, variableIndex, floorValue, true, $"? {floorValue}"));   // x ? floor
                children.Add(CreateChildWithConstraint(parentNode, variableIndex, ceilValue, false, $"? {ceilValue}"));    // x ? ceil
            }

            return children;
        }

        /// <summary>
        /// Creates a child node with an additional constraint
        /// </summary>
        private BranchAndBoundNode CreateChildWithConstraint(BranchAndBoundNode parent, int variableIndex, double boundValue, bool isUpperBound, string constraintDescription)
        {
            var childModel = CloneModel(parent.Model);
            
            // Add the branching constraint to the model
            AddBranchingConstraint(childModel, variableIndex, boundValue, isUpperBound);
            
            return new BranchAndBoundNode
            {
                Id = nextNodeId++,
                ParentId = parent.Id,
                Model = childModel,
                BranchingConstraint = $"x{variableIndex + 1} {constraintDescription}",
                BranchedVariableIndex = variableIndex,
                Depth = parent.Depth + 1
            };
        }

        /// <summary>
        /// Adds a branching constraint to the canonical model
        /// </summary>
        private void AddBranchingConstraint(CanonicalLinearProgrammingModel model, int variableIndex, double boundValue, bool isUpperBound)
        {
            int numVariables = model.CoefficientMatrix[0].Length;
            int numConstraints = model.CoefficientMatrix.Length;

            // Create new constraint coefficient array
            double[] newConstraint = new double[numVariables];
            newConstraint[variableIndex] = isUpperBound ? 1 : -1; // For x ? b: +1x; For x ? b: -1x

            // Create new RHS value
            double newRHS = isUpperBound ? boundValue : -boundValue; // For x ? b: b; For x ? b: -b

            // Expand the coefficient matrix
            var newMatrix = new double[numConstraints + 1][];
            for (int i = 0; i < numConstraints; i++)
            {
                newMatrix[i] = (double[])model.CoefficientMatrix[i].Clone();
            }
            newMatrix[numConstraints] = newConstraint;

            // Expand the RHS vector
            var newRHSVector = new double[numConstraints + 1];
            for (int i = 0; i < numConstraints; i++)
            {
                newRHSVector[i] = model.RHSVector[i];
            }
            newRHSVector[numConstraints] = newRHS;

            // Update the model
            model.CoefficientMatrix = newMatrix;
            model.RHSVector = newRHSVector;
        }

        /// <summary>
        /// Creates a deep copy of a canonical model
        /// </summary>
        private CanonicalLinearProgrammingModel CloneModel(CanonicalLinearProgrammingModel original)
        {
            return new CanonicalLinearProgrammingModel
            {
                CoefficientMatrix = original.CoefficientMatrix.Select(row => (double[])row.Clone()).ToArray(),
                RHSVector = (double[])original.RHSVector.Clone(),
                ObjectiveCoefficients = (double[])original.ObjectiveCoefficients.Clone(),
                VariableTypes = (VariableType[])original.VariableTypes.Clone()
            };
        }

        /// <summary>
        /// Utility method to print solution vector
        /// </summary>
        private void PrintSolution(double[] solution)
        {
            for (int i = 0; i < solution.Length; i++)
            {
                Console.WriteLine($"x{i + 1} = {solution[i]:F3}");
            }
        }

        /// <summary>
        /// Displays the canonical form for debugging
        /// </summary>
        public void DisplayCanonicalForm(CanonicalLinearProgrammingModel model)
        {
            Console.WriteLine("=== CANONICAL FORM ===");
            
            // Display objective function
            Console.Write("max  ");
            for (int j = 0; j < model.ObjectiveCoefficients.Length; j++)
            {
                Console.Write($"{model.ObjectiveCoefficients[j]:+0.000;-0.000}x{j + 1} ");
            }
            Console.WriteLine();
            
            // Display constraints
            Console.WriteLine("subject to:");
            for (int i = 0; i < model.CoefficientMatrix.Length; i++)
            {
                Console.Write("     ");
                for (int j = 0; j < model.CoefficientMatrix[i].Length; j++)
                {
                    Console.Write($"{model.CoefficientMatrix[i][j]:+0.000;-0.000}x{j + 1} ");
                }
                Console.WriteLine($"= {model.RHSVector[i]:F3}");
            }
            
            // Display variable bounds
            Console.WriteLine("where:");
            for (int j = 0; j < model.VariableTypes.Length; j++)
            {
                string bounds = model.VariableTypes[j] switch
                {
                    VariableType.NonNegative => "? 0",
                    VariableType.NonPositive => "? 0", 
                    VariableType.Unrestricted => "unrestricted",
                    VariableType.Integer => "? 0, integer",
                    VariableType.Binary => "? {0,1}",
                    _ => "unknown"
                };
                Console.WriteLine($"     x{j + 1} {bounds}");
            }
            Console.WriteLine();
        }
    }
}