namespace LinearProgramming.Parsing
{
    /// <summary>
    /// Specifies the type of constraint in a linear programming problem.
    /// </summary>
    public enum ConstraintType
    {
        /// <summary>
        /// Less than or equal to (≤) constraint
        /// </summary>
        LessThanOrEqual,
        
        /// <summary>
        /// Equal to (=) constraint
        /// </summary>
        Equal,
        
        /// <summary>
        /// Greater than or equal to (≥) constraint
        /// </summary>
        GreaterThanOrEqual
    }
}
