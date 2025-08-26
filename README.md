# Linear Programming and Integer Programming Solver

A comprehensive .NET 9 application that solves Linear Programming (LP) and Integer Programming (IP) models with sensitivity analysis capabilities.

## Project Overview

This Visual Studio project implements various optimization algorithms for solving linear programming and integer programming problems. The application accepts mathematical models via input text files and exports detailed results to output files.

## Architecture

The solution consists of four main projects:

- **LinearProgrammingApp** - Console application and main entry point
- **LinearProgramming.Parsing** - Input file parsing and model representation
- **LinearProgramming.Algorithms** - Implementation of optimization algorithms
- **LinearProgramming.Tests** - Unit tests for validation

## ? Completed Features

### Input/Output System
- ? **Universal Input Parser** - Accepts text files with mathematical models in the specified format
- ? **Flexible Variable Support** - Handles random amount of decision variables
- ? **Flexible Constraint Support** - Handles random amount of constraints
- ? **Variable Type Support** - Supports +, -, urs, int, bin variable types
- ? **Constraint Type Support** - Handles <=, >=, = constraints
- ? **Error Handling** - Comprehensive input validation and error reporting

### Core Algorithms
- ? **Primal Simplex Algorithm** - Full implementation with tableau-based approach
- ? **Revised Primal Simplex Algorithm** - Matrix-based implementation with Phase I/II method
- ? **Canonical Form Conversion** - Automatic conversion from input format to canonical form
- ? **Artificial Variables Handling** - Proper management of artificial variables for = and >= constraints
- ? **Special Cases Detection** - Identifies and handles infeasible and unbounded solutions

### Programming Quality
- ? **Comprehensive Comments** - Well-documented code with XML documentation
- ? **Best Practices Implementation** - SOLID principles, proper error handling, and clean architecture
- ? **Object-Oriented Design** - Modular, extensible design pattern
- ? **Unit Testing Framework** - XUnit tests for validation

### Mathematical Correctness
- ? **Phase I/II Method** - Proper handling of artificial variables
- ? **Degeneracy Handling** - Robust handling of degenerate solutions
- ? **Numerical Stability** - Epsilon tolerance for floating-point operations
- ? **Basis Management** - Intelligent basis identification and repair

### Output Capabilities
- ? **Solution Display** - Shows optimal solution vector and objective value
- ? **Status Reporting** - Reports Optimal, Infeasible, Unbounded, or iteration limit status
- ? **Tableau History** - Stores all simplex tableau iterations

## Input File Format

The application accepts input files in the following format:

```
max +2 +3 +3 +5 +2 +4
+11 +8 +6 +14 +10 +10 <= +40
bin bin bin bin bin bin
```

**Format Specification:**
- Line 1: `max/min` followed by objective coefficients with +/- operators
- Lines 2-n: Constraint coefficients, operator (<=, >=, =), and RHS value
- Last line: Variable types (+, -, urs, int, bin) for each variable

## Usage

```bash
LinearProgrammingApp.exe <inputfile>
```

The program will prompt you to choose between:
1. Primal Simplex Algorithm
2. Revised Primal Simplex Algorithm

## ? Still To Be Implemented

### Integer Programming Algorithms
- ? **Branch and Bound Simplex Algorithm** - For integer programming problems
- ? **Cutting Plane Algorithm** - Alternative IP solving method  
- ? **Branch and Bound Knapsack Algorithm** - Specialized for knapsack problems

### Advanced Features
- ? **Menu-Driven Interface** - Interactive console menu system
- ? **Output File Export** - Automatic export of results to text files
- ? **Tableau Display** - Show all simplex tableau iterations in output
- ? **Product Form Display** - Show B?¹ matrix updates for revised simplex
- ? **Decimal Precision Control** - Round all values to 3 decimal places

### Sensitivity Analysis
- ? **Non-Basic Variable Range Analysis** - Display and apply changes
- ? **Basic Variable Range Analysis** - Display and apply changes  
- ? **RHS Range Analysis** - Constraint right-hand-side sensitivity
- ? **Coefficient Change Analysis** - Technology coefficient sensitivity
- ? **New Activity Addition** - Add new variables to optimal solution
- ? **New Constraint Addition** - Add new constraints to optimal solution
- ? **Shadow Prices Calculation** - Display dual variable values

### Duality Analysis  
- ? **Dual Problem Generation** - Create and solve the dual formulation
- ? **Strong/Weak Duality Verification** - Mathematical duality validation
- ? **Complementary Slackness** - Verify optimality conditions

### Enhanced Algorithm Features
- ? **Backtracking Implementation** - For branch and bound algorithms
- ? **Sub-problem Generation** - Create all possible branching nodes
- ? **Node Fathoming** - Eliminate sub-optimal branches
- ? **Best Candidate Tracking** - Monitor best solution across branches

## Current Status

**Foundation Complete (60%)** - The core linear programming solving capabilities are fully implemented and mathematically correct. The application can successfully solve LP problems using both standard and revised simplex methods.

**Next Phase** - Focus on implementing integer programming algorithms, sensitivity analysis features, and enhanced user interface components.

## Technical Requirements

- **.NET 9.0** - Latest .NET framework
- **Visual Studio 2022** - Recommended IDE
- **MathNet.Numerics** - Mathematical operations library
- **XUnit** - Testing framework

## Building and Running

```bash
# Build the solution
dotnet build

# Run with input file
dotnet run --project LinearProgrammingApp -- input.txt
```

## Testing

```bash
# Run all tests
dotnet test
```

---

**Note**: This is an academic project for LPR 381 (Linear Programming) course. The implementation prioritizes mathematical correctness and educational value over performance optimization.
