namespace LinearProgramming.Algorithms
{
    public class Variable
    {
        public string Name { get; set; } = string.Empty;
        public double Coefficient { get; set; }
        public double LowerBound { get; set; } = 0;
        public double UpperBound { get; set; } = double.PositiveInfinity;
        public bool IsInteger { get; set; }
    }
}
