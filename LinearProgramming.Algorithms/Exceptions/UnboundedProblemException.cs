using System;

namespace LinearProgramming.Algorithms.Exceptions
{
    public class UnboundedProblemException : Exception
    {
        public UnboundedProblemException() : base("The linear program is unbounded.")
        {
        }

        public UnboundedProblemException(string message) : base(message)
        {
        }

        public UnboundedProblemException(string message, Exception innerException) 
            : base(message, innerException)
        {
        }
    }
}
