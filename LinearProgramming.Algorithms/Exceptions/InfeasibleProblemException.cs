using System;

namespace LinearProgramming.Algorithms.Exceptions
{
    public class InfeasibleProblemException : Exception
    {
        public InfeasibleProblemException() : base("The linear program is infeasible.")
        {
        }

        public InfeasibleProblemException(string message) : base(message)
        {
        }

        public InfeasibleProblemException(string message, Exception innerException) 
            : base(message, innerException)
        {
        }
    }
}
