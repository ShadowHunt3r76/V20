using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace LPRProject
{
    public class Knapsack
    {
        public double SolveKnapsack(double[] weights, double[] values, double capacity)
        {
            int n = weights.Length;
            var dp = Matrix<double>.Build.Dense(n + 1, (int)capacity + 1);

            for (int i = 0; i <= n; i++)
            {
                for (int w = 0; w <= capacity; w++)
                {
                    if (i == 0 || w == 0)
                    {
                        dp[i, w] = 0;
                    }
                    else if (weights[i - 1] <= w)
                    {
                        dp[i, w] = Math.Max(values[i - 1] + dp[i - 1, w - (int)weights[i - 1]], dp[i - 1, w]);
                    }
                    else
                    {
                        dp[i, w] = dp[i - 1, w];
                    }
                }
            }

            return dp[n, (int)capacity];
        }

        public static void Main(string[] args)
        {
            try
            {
                Console.WriteLine("Enter the number of items:");
                int n = int.Parse(Console.ReadLine());

                if (n <= 0)
                {
                    throw new ArgumentException("Number of items must be positive.");
                }

                double[] weights = new double[n];
                double[] values = new double[n];

                for (int i = 0; i < n; i++)
                {
                    Console.WriteLine($"Enter weight of item {i + 1}:");
                    weights[i] = double.Parse(Console.ReadLine());

                    if (weights[i] < 0)
                    {
                        throw new ArgumentException("Weight must be non-negative.");
                    }

                    Console.WriteLine($"Enter value of item {i + 1}:");
                    values[i] = double.Parse(Console.ReadLine());

                    if (values[i] < 0)
                    {
                        throw new ArgumentException("Value must be non-negative.");
                    }
                }

                Console.WriteLine("Enter the capacity of the knapsack:");
                double capacity = double.Parse(Console.ReadLine());

                if (capacity <= 0)
                {
                    throw new ArgumentException("Capacity must be positive.");
                }

                Knapsack knapsackSolver = new Knapsack();
                double maxValue = knapsackSolver.SolveKnapsack(weights, values, capacity);

                Console.WriteLine($"Maximum value in knapsack = {maxValue}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
            }
        }
    }
}