using System;
using System.Collections.Generic;
using System.Linq;
using LinearProgramming.Parsing;
using LinearProgramming.Algorithms.PrimalSimplex;
using LinearProgramming.Algorithms.Utils;
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

            // Display the complete branch and bound tree
            DisplayBranchAndBoundTree(solution);

            // Display the path to the best solution
            if (solution.BestCandidate != null)
            {
                DisplaySolutionPath(solution);
            }
            
            var results = new List<object[]>
            {
                new object[] { "Total Nodes Explored", solution.TotalNodesExplored },
                new object[] { "Total Nodes Fathomed", solution.TotalNodesFathomed },
                new object[] { "Best Integer Solution", solution.BestObjectiveValue.ToString("F3") },
                new object[] { "Status", solution.OverallStatus }
            };

            var table = OutputFormatter.CreateTable(
                new[] { "Metric", "Value" },
                results
            );

            Console.WriteLine(OutputFormatter.CreateBox(table, "BRANCH AND BOUND RESULTS"));

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
        /// <summary>
        /// Analyzes infeasibility and provides detailed reason
        /// </summary>
        private string AnalyzeInfeasibility(CanonicalLinearProgrammingModel model, LinearProgramSolution solution)
        {
            var analysis = new StringBuilder();
            
            // Check for trivially infeasible constraints
            for (int i = 0; i < model.ConstraintTypes.Length; i++)
            {
                bool allNonPositive = true;
                bool allNonNegative = true;
                bool allZero = true;
                bool hasPositiveRHS = model.RHSVector[i] > EPSILON;
                bool hasNegativeRHS = model.RHSVector[i] < -EPSILON;

                foreach (var coef in model.CoefficientMatrix[i])
                {
                    if (coef > EPSILON)
                    {
                        allNonPositive = false;
                        allZero = false;
                    }
                    else if (coef < -EPSILON)
                    {
                        allNonNegative = false;
                        allZero = false;
                    }
                }

                // Check for trivially infeasible constraints
                if ((model.ConstraintTypes[i] == ConstraintType.GreaterThanOrEqual && allNonPositive && hasPositiveRHS) ||
                    (model.ConstraintTypes[i] == ConstraintType.LessThanOrEqual && allNonNegative && hasNegativeRHS) ||
                    (model.ConstraintTypes[i] == ConstraintType.Equal && allZero && Math.Abs(model.RHSVector[i]) > EPSILON))
                {
                    analysis.AppendLine($"- Constraint {i + 1} is trivially infeasible");
                    analysis.AppendLine($"  {FormatConstraint(model.CoefficientMatrix[i], model.ConstraintTypes[i], model.RHSVector[i], model.VariableNames)}");
                }
            }

            // Check for conflicting constraints
            for (int i = 0; i < model.ConstraintTypes.Length; i++)
            {
                for (int j = i + 1; j < model.ConstraintTypes.Length; j++)
                {
                    if (AreConflictingConstraints(
                        model.CoefficientMatrix[i], model.ConstraintTypes[i], model.RHSVector[i],
                        model.CoefficientMatrix[j], model.ConstraintTypes[j], model.RHSVector[j],
                        out string conflictReason))
                    {
                        analysis.AppendLine($"- Conflict between constraints {i + 1} and {j + 1}: {conflictReason}");
                    }
                }
            }

            return analysis.Length > 0 ? analysis.ToString() : "No specific infeasibility reason identified";
        }

        private bool AreConflictingConstraints(
            double[] coeffs1, ConstraintType type1, double rhs1,
            double[] coeffs2, ConstraintType type2, double rhs2,
            out string reason)
        {
            // Check if constraints are parallel but have no feasible region
            bool allEqual = true;
            for (int k = 0; k < coeffs1.Length; k++)
            {
                if (Math.Abs(coeffs1[k] - coeffs2[k]) > EPSILON)
                {
                    allEqual = false;
                    break;
                }
            }
            
            if (allEqual)
            {
                if ((type1 == ConstraintType.LessThanOrEqual && type2 == ConstraintType.GreaterThanOrEqual && rhs1 < rhs2) ||
                    (type1 == ConstraintType.GreaterThanOrEqual && type2 == ConstraintType.LessThanOrEqual && rhs1 > rhs2))
                {
                    reason = "Parallel constraints with no feasible region";
                    return true;
                }
            }
            
            reason = string.Empty;
            return false;
        }

        private bool ShouldFathomNode(BranchAndBoundNode node, double bestKnownValue, ParsedLinearProgrammingModel originalModel)
        {
            // Fathom by infeasibility
            if (node.Solution.Status == "Infeasible")
            {
                node.FathomReason = "Infeasible: " + AnalyzeInfeasibility(node.Model, node.Solution);
                Console.WriteLine($"\n=== Infeasibility Analysis for Node {node.Id} ===");
                Console.WriteLine(node.FathomReason);
                return true;
            }

            // Fathom by unboundedness
            if (node.Solution.Status == "Unbounded")
            {
                node.FathomReason = "Unbounded LP relaxation - no optimal integer solution in this branch";
                Console.WriteLine($"\n=== Unbounded Analysis for Node {node.Id} ===");
                Console.WriteLine("The LP relaxation is unbounded. This suggests that the integer problem");
                Console.WriteLine("may be unbounded or may require additional constraints to find an optimal solution.");
                return true;
            }

            // Fathom by bound (for maximization: if upper bound <= best known integer value)
            if (node.Solution.Status == "Optimal")
            {
                node.UpperBound = node.Solution.ObjectiveValue;
                
                // Check for numerical stability issues
                if (double.IsInfinity(node.UpperBound) || double.IsNaN(node.UpperBound))
                {
                    node.FathomReason = "Numerical instability detected in solution";
                    return true;
                }
                
                // Check if this node's bound is worse than the best known solution
                if (bestKnownValue != double.MinValue && node.UpperBound <= bestKnownValue + EPSILON)
                {
                    node.FathomReason = $"Fathomed by bound: {node.UpperBound:F3} (current best: {bestKnownValue:F3})";
                    Console.WriteLine($"\n=== Bound Analysis for Node {node.Id} ===");
                    Console.WriteLine($"Node {node.Id} has upper bound {node.UpperBound:F3} which is not better than current best {bestKnownValue:F3}");
                    return true;
                }
                
                // Check for near-integer solution that's not quite feasible
                if (IsNearlyIntegerFeasible(node.Solution.SolutionVector, originalModel, out var fractionalVar, out var fraction))
                {
                    Console.WriteLine($"\n=== Near-Integer Solution in Node {node.Id} ===");
                    Console.WriteLine($"Variable x{fractionalVar + 1} has value {node.Solution.SolutionVector[fractionalVar]:F6} (fractional part: {fraction:F6})");
                    Console.WriteLine("Consider adjusting the branching strategy or tolerance levels.");
                }
            }

            return false;
        }

        /// <summary>
        /// Checks if the solution satisfies integer constraints
        /// </summary>
        /// <summary>
        /// Checks if a solution is integer feasible within tolerance
        /// </summary>
        private bool IsIntegerFeasible(double[] solution, ParsedLinearProgrammingModel originalModel)
        {
            if (solution == null) 
                return false;

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
        /// Checks if a solution is nearly integer feasible and returns the most fractional variable
        /// </summary>
        private bool IsNearlyIntegerFeasible(double[] solution, ParsedLinearProgrammingModel originalModel, out int fractionalVar, out double fraction)
        {
            fractionalVar = -1;
            fraction = 0.0;
            
            if (solution == null) 
                return false;
                
            double maxFraction = 0.0;
            
            for (int i = 0; i < Math.Min(solution.Length, originalModel.Variables.Count); i++)
            {
                var varType = originalModel.Variables[i].Type;
                if (varType == VariableType.Continuous)
                    continue;
                    
                double value = solution[i];
                double intPart = Math.Floor(value);
                double frac = Math.Abs(value - Math.Round(value));
                
                // Consider values very close to integers as integral
                if (frac > EPSILON && frac < 1.0 - EPSILON)
                {
                    if (frac > maxFraction)
                    {
                        maxFraction = frac;
                        fractionalVar = i;
                        fraction = frac;
                    }
                }
            }
            
            return fractionalVar >= 0;
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
        /// Performs sensitivity analysis on the solution
        /// </summary>
        /// <param name="solution">The branch and bound solution to analyze</param>
        /// <param name="originalModel">The original parsed model</param>
        public static void PerformSensitivityAnalysis(BranchAndBoundSolution solution, ParsedLinearProgrammingModel originalModel)
        {
            if (solution.BestCandidate != null && solution.BestCandidate.Solution?.Status == "Optimal")
            {
                try
                {
                    var sensitivityAnalysis = new BranchAndBoundSensitivityAnalysis(solution, originalModel);
                    sensitivityAnalysis.PerformAnalysis();
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"\n⚠ Warning: Could not perform sensitivity analysis: {ex.Message}");
                    Console.WriteLine("The solution is still valid, but sensitivity analysis is not available.");
                }
            }
            else
            {
                Console.WriteLine("\n⚠ Cannot perform sensitivity analysis: No optimal solution found.");
            }
        }

        /// <summary>
        /// Displays the branch and bound tree structure
        /// </summary>
        private void DisplayBranchAndBoundTree(BranchAndBoundSolution solution)
        {
            Console.WriteLine("\n╔══════════════════════════════════════════════════════╗");
            Console.WriteLine("║               BRANCH AND BOUND TREE                   ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════╝");
            
            // Group nodes by depth
            var nodesByDepth = solution.AllNodes.GroupBy(n => n.Depth).OrderBy(g => g.Key);
            
            foreach (var depthGroup in nodesByDepth)
            {
                Console.Write($"\nDepth {depthGroup.Key}: ");
                
                foreach (var node in depthGroup.OrderBy(n => n.Id))
                {
                    Console.ForegroundColor = GetNodeColor(node, solution.BestCandidate);
                    string nodeStatus = GetNodeStatusSymbol(node);
                    Console.Write($"[{nodeStatus}{node.Id}] ");
                    Console.ResetColor();
                }
                
                // Draw connecting lines for the next level
                if (depthGroup.Key < nodesByDepth.Max(g => g.Key))
                {
                    Console.WriteLine("\n" + new string(' ', 8) + 
                        string.Join("   ", 
                            depthGroup.SelectMany(n => 
                                solution.AllNodes.Where(c => c.ParentId == n.Id)
                                    .Select(c => $"│  ")
                            )
                        )
                    );
                    
                    Console.WriteLine(new string(' ', 8) + 
                        string.Join("   ", 
                            depthGroup.SelectMany(n => 
                                solution.AllNodes.Where(c => c.ParentId == n.Id)
                                    .Select(c => n.Id == solution.AllNodes.Last(n2 => n2.Depth == depthGroup.Key).Id ? "└─" : "├─")
                            )
                        )
                    );
                }
            }
            
            // Display legend
            Console.WriteLine("\n\n╔══════════════════════════════════════════════════════╗");
            Console.WriteLine("║                     LEGEND                            ║");
            Console.WriteLine("╠══════════════════════════════════════════════════════╣");
            Console.Write("║ ");
            Console.ForegroundColor = ConsoleColor.Green;
            Console.Write("●");
            Console.ResetColor();
            Console.WriteLine(" - Best Integer Solution  ");
            Console.Write("║ ");
            Console.ForegroundColor = ConsoleColor.Blue;
            Console.Write("○");
            Console.ResetColor();
            Console.WriteLine(" - Active Node  ");
            Console.Write("║ ");
            Console.ForegroundColor = ConsoleColor.Red;
            Console.Write("✗");
            Console.ResetColor();
            Console.WriteLine(" - Fathomed Node  ");
            Console.Write("║ ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.Write("◌");
            Console.ResetColor();
            Console.WriteLine(" - Unexplored Node  ");
            Console.WriteLine("╚══════════════════════════════════════════════════════╝");
        }
        
        /// <summary>
        /// Displays the path from root to the best solution
        /// </summary>
        private void DisplaySolutionPath(BranchAndBoundSolution solution)
        {
            Console.WriteLine("\n╔══════════════════════════════════════════════════════╗");
            Console.WriteLine("║               PATH TO BEST SOLUTION                   ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════╝");
            
            // Reconstruct the path from best node to root
            var path = new List<BranchAndBoundNode>();
            var currentNode = solution.BestCandidate;
            
            while (currentNode != null)
            {
                path.Insert(0, currentNode);
                currentNode = solution.AllNodes.FirstOrDefault(n => n.Id == currentNode.ParentId);
            }
            
            // Display the path
            for (int i = 0; i < path.Count; i++)
            {
                var node = path[i];
                string prefix = new string(' ', i * 2);
                string connector = i > 0 ? "└─ " : "";
                
                Console.ForegroundColor = GetNodeColor(node, solution.BestCandidate);
                Console.Write($"{prefix}{connector}[{node.Id}] ");
                Console.ResetColor();
                
                if (i > 0)
                {
                    Console.WriteLine($"{node.BranchingConstraint}");
                }
                else
                {
                    Console.WriteLine("Root");
                }
                
                if (i == path.Count - 1) // Last node (best solution)
                {
                    Console.WriteLine($"{new string(' ', (i + 1) * 2 + 2)}Objective: {node.Solution.ObjectiveValue:F3}");
                    Console.WriteLine($"{new string(' ', (i + 1) * 2 + 2)}Variables: {string.Join(", ", node.Solution.SolutionVector.Select((v, idx) => $"x{idx + 1} = {v:F3}"))}");
                }
            }
        }
        
        /// <summary>
        /// Gets the appropriate color for a node based on its status
        /// </summary>
        private ConsoleColor GetNodeColor(BranchAndBoundNode node, BranchAndBoundNode bestNode)
        {
            if (node == bestNode) return ConsoleColor.Green;
            if (node.IsFathomed) return ConsoleColor.Red;
            if (node.Solution != null && node.Solution.Status == "Optimal") return ConsoleColor.Blue;
            return ConsoleColor.Gray;
        }
        
        /// <summary>
        /// Gets a symbol representing the node's status
        /// </summary>
        private string GetNodeStatusSymbol(BranchAndBoundNode node)
        {
            if (node.IsFathomed) return "✗";
            if (node.Solution != null && node.Solution.Status == "Optimal") return "●";
            return "◌";
        }
        
        /// <summary>
        /// Utility method to print solution vector with formatting
        /// </summary>
        private void PrintSolution(double[] solution)
        {
            if (solution == null)
            {
                Console.WriteLine("No solution found.");
                return;
            }

            Console.WriteLine("\n╔══════════════════════════════════════════════════════╗");
            Console.WriteLine("║                   SOLUTION VECTOR                    ║");
            Console.WriteLine("╠══════════════════════════════════════════════════════╣");
            
            for (int i = 0; i < solution.Length; i++)
            {
                Console.WriteLine($"║ x{i + 1,-4} = {solution[i],-10:F6} {new string(' ', 30 - solution[i].ToString("F6").Length)}║");
            }
            
            Console.WriteLine("╚══════════════════════════════════════════════════════╝");
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