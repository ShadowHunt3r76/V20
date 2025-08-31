namespace LinearProgramming.Parsing
{
    /// <summary>
    /// Specifies the type and bounds for variables in a linear programming model.
    /// </summary>
    public enum VariableType
    {
        /// <summary>
        /// Non-negative variable (≥ 0)
        /// </summary>
        NonNegative,
        
        /// <summary>
        /// Non-positive variable (≤ 0)
        /// </summary>
        NonPositive,
        
        /// <summary>
        /// Unrestricted variable (can be positive, negative, or zero)
        /// </summary>
        Unrestricted,
        
        /// <summary>
        /// Integer variable (whole numbers only)
        /// </summary>
        Integer,
        
        /// <summary>
        /// Binary variable (0 or 1 only)
        /// </summary>
        Binary
    }
}
