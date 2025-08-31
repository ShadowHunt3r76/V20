using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LinearProgramming.Algorithms.PrimalSimplex
{
    /// <summary>
    /// Provides common matrix and vector operations used by simplex solvers
    /// </summary>
    public static class MatrixUtils
    {
        public const double Epsilon = 1e-10;

        /// <summary>
        /// Computes the dot product of two vectors
        /// </summary>
        public static double DotProduct(double[] a, double[] b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (a.Length != b.Length) 
                throw new ArgumentException("Vectors must have the same length");
                
            double result = 0;
            for (int i = 0; i < a.Length; i++)
            {
                result += a[i] * b[i];
            }
            return result;
        }

        /// <summary>
        /// Creates an identity matrix of the specified size
        /// </summary>
        public static double[,] CreateIdentityMatrix(int size)
        {
            if (size <= 0) throw new ArgumentOutOfRangeException(nameof(size), "Size must be positive");
            
            var identity = new double[size, size];
            for (int i = 0; i < size; i++)
            {
                identity[i, i] = 1.0;
            }
            return identity;
        }

        /// <summary>
        /// Converts a 1D array to a 2D array with a single row
        /// </summary>
        public static double[,] To2D(double[] array)
        {
            if (array == null) throw new ArgumentNullException(nameof(array));
            
            var result = new double[1, array.Length];
            for (int i = 0; i < array.Length; i++)
            {
                result[0, i] = array[i];
            }
            return result;
        }

        /// <summary>
        /// Converts a jagged array to a 2D array
        /// </summary>
        /// <param name="jaggedArray">The jagged array to convert</param>
        /// <returns>A 2D array with the same elements as the input jagged array</returns>
        public static double[,] ConvertToMatrix(double[][] jaggedArray)
        {
            if (jaggedArray == null || jaggedArray.Length == 0)
                return new double[0, 0];

            int rows = jaggedArray.Length;
            int cols = jaggedArray[0].Length;
            double[,] matrix = new double[rows, cols];
            
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    matrix[i, j] = jaggedArray[i][j];
                }
            }
            return matrix;
        }


        public static double[,] InvertMatrix(double[,] matrix, double epsilon = 1e-10)
        {
            int n = matrix.GetLength(0);
            double[,] inverse = new double[n, n];
            
            // Initialize inverse as identity matrix
            for (int i = 0; i < n; i++)
                inverse[i, i] = 1.0;
                
            // Create a copy of the matrix to avoid modifying the original
            double[,] tempMatrix = (double[,])matrix.Clone();
                
            // Perform Gaussian elimination
            for (int col = 0; col < n; col++)
            {
                // Find pivot row
                int pivot = col;
                double maxVal = Math.Abs(tempMatrix[col, col]);
                
                for (int row = col + 1; row < n; row++)
                {
                    if (Math.Abs(tempMatrix[row, col]) > maxVal)
                    {
                        maxVal = Math.Abs(tempMatrix[row, col]);
                        pivot = row;
                    }
                }
                
                // Check for singular matrix
                if (maxVal < epsilon)
                    return null;
                    
                // Swap rows if needed
                if (pivot != col)
                {
                    for (int i = 0; i < n; i++)
                    {
                        // Swap in original matrix
                        double temp = tempMatrix[col, i];
                        tempMatrix[col, i] = tempMatrix[pivot, i];
                        tempMatrix[pivot, i] = temp;
                        
                        // Swap in inverse matrix
                        temp = inverse[col, i];
                        inverse[col, i] = inverse[pivot, i];
                        inverse[pivot, i] = temp;
                    }
                }
                
                // Normalize pivot row
                double pivotVal = tempMatrix[col, col];
                for (int i = 0; i < n; i++)
                {
                    tempMatrix[col, i] /= pivotVal;
                    inverse[col, i] /= pivotVal;
                }
                
                // Eliminate other rows
                for (int row = 0; row < n; row++)
                {
                    if (row != col && Math.Abs(tempMatrix[row, col]) > epsilon)
                    {
                        double factor = tempMatrix[row, col];
                        for (int i = 0; i < n; i++)
                        {
                            tempMatrix[row, i] -= factor * tempMatrix[col, i];
                            inverse[row, i] -= factor * inverse[col, i];
                        }
                    }
                }
            }
            
            return inverse;
        }
        
        public static double[] MatrixVectorMultiply(double[,] matrix, double[] vector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[] result = new double[rows];
            
            for (int i = 0; i < rows; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < cols; j++)
                {
                    sum += matrix[i, j] * vector[j];
                }
                result[i] = sum;
            }
            
            return result;
        }
        
        public static double[] VectorMatrixMultiply(double[] vector, double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[] result = new double[cols];
            
            for (int j = 0; j < cols; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < rows; i++)
                {
                    sum += vector[i] * matrix[i, j];
                }
                result[j] = sum;
            }
            
            return result;
        }
        
        public static double[,] ExtractBasisMatrix(double[][] A, int[] basisIndices, int size)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (basisIndices == null) throw new ArgumentNullException(nameof(basisIndices));
            
            double[,] basisMatrix = new double[size, size];
            
            for (int i = 0; i < size; i++)
            {
                int col = basisIndices[i];
                for (int j = 0; j < size; j++)
                {
                    basisMatrix[j, i] = A[j][col];
                }
            }
            
            return basisMatrix;
        }
        
        public static double[] GetColumn(double[][] matrix, int colIndex)
        {
            if (matrix == null) throw new ArgumentNullException(nameof(matrix));
            if (colIndex < 0 || colIndex >= matrix[0].Length)
                throw new ArgumentOutOfRangeException(nameof(colIndex));
                
            double[] column = new double[matrix.Length];
            for (int i = 0; i < matrix.Length; i++)
            {
                column[i] = matrix[i][colIndex];
            }
            return column;
        }
        
        public static T[,] To2D<T>(T[][] source)
        {
            if (source == null) throw new ArgumentNullException(nameof(source));
            
            int rows = source.Length;
            if (rows == 0) return new T[0, 0];
            
            int cols = source[0].Length;
            T[,] result = new T[rows, cols];
            
            for (int i = 0; i < rows; i++)
            {
                if (source[i] == null || source[i].Length != cols)
                    throw new ArgumentException("Jagged array is not rectangular", nameof(source));
                    
                for (int j = 0; j < cols; j++)
                {
                    result[i, j] = source[i][j];
                }
            }
            
            return result;
        }
        /// <summary>
        /// Formats a tableau for display
        /// </summary>
        public static string FormatTableau(double[,] tableau)
        {
            if (tableau == null) return "[null tableau]";
            
            var sb = new StringBuilder();
            int rows = tableau.GetLength(0);
            int cols = tableau.GetLength(1);
            
            // Determine column widths
            var colWidths = new int[cols];
            for (int j = 0; j < cols; j++)
            {
                int maxLen = 0;
                for (int i = 0; i < rows; i++)
                {
                    string s = string.Format("{0,10:0.#####}", tableau[i, j]);
                    maxLen = Math.Max(maxLen, s.Length);
                }
                colWidths[j] = maxLen + 1;
            }
            
            // Print the tableau
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    sb.Append($"{tableau[i, j],10:0.#####} ");
                }
                sb.AppendLine();
            }
            
            return sb.ToString();
        }
    }
}
