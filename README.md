# Linear Programming and Integer Programming Solver

A .NET 9 application that solves Linear Programming (LP) and basic Integer Programming (IP) models. This project was developed as part of the LPR 381 (Linear Programming) course, focusing on implementing core optimization algorithms with an emphasis on educational value and code clarity.

## Project Overview

This Visual Studio project implements various optimization algorithms for solving linear programming and integer programming problems. The application accepts mathematical models via input text files and exports detailed results to output files.

## Architecture

The solution consists of four main projects:

- **LinearProgrammingApp** - Console application and main entry point
- **LinearProgramming.Parsing** - Input file parsing and model representation
- **LinearProgramming.Algorithms** - Implementation of optimization algorithms
- **LinearProgramming.Tests** - Unit tests for validation

## ✅ Implemented Features

### Input/Output System
- **Universal Input Parser** - Accepts text files with mathematical models in the specified format
- **Flexible Variable Support** - Handles any number of decision variables
- **Flexible Constraint Support** - Handles any number of constraints
- **Variable Type Support** - Supports +, -, urs, int, bin variable types
- **Constraint Type Support** - Handles <=, >=, = constraints
- **Error Handling** - Input validation and error reporting

### Core Algorithms
- **Primal Simplex Algorithm** - Tableau-based implementation
- **Revised Primal Simplex Algorithm** - Matrix-based implementation with Phase I/II method
- **Canonical Form Conversion** - Automatic conversion from input format to canonical form
- **Artificial Variables** - Management of artificial variables for = and >= constraints
- **Special Cases** - Handles infeasible and unbounded solutions

### Code Quality
- **Documentation** - XML documentation for public APIs
- **Error Handling** - Basic error handling and input validation
- **Object-Oriented Design** - Follows basic OOP principles
- **Unit Tests** - Basic test coverage for core functionality

### Mathematical Implementation
- **Phase I/II Method** - Basic implementation for artificial variables
- **Numerical Stability** - Basic floating-point tolerance handling
- **Basis Management** - Basic basis handling

### Output
- **Solution Display** - Shows solution vector and objective value
- **Status Reporting** - Basic status reporting (Optimal/Infeasible/Unbounded)

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

## 🚧 In Progress / To Be Implemented

### Branch and Bound Implementation (In Progress)
- **Basic Structure** - Initial implementation started
- **Node Management** - Basic node creation and management
- **Branching** - Basic variable selection and branching
- **Bound Calculation** - Basic bound calculation
- **Pruning** - Basic pruning of suboptimal branches

### Future Enhancements
- **Cutting Plane Algorithm** - For integer programming problems
- **Knapsack Solver** - Specialized branch and bound for knapsack problems

### Potential Future Features
- **Enhanced Output** - Better formatting and display of solutions
- **File I/O** - Export/import of problems and solutions
- **Precision Control** - Configurable decimal precision
- **Performance Optimizations** - For larger problem instances

### Advanced Topics (Future Work)
- **Sensitivity Analysis** - RHS and coefficient sensitivity
- **Duality** - Dual problem generation and analysis
- **Algorithm Enhancements** - Advanced branching rules and heuristics

## Current Status

**Core Functionality (80%)** - The application can solve standard linear programming problems using both primal and revised simplex methods. The basic structure for branch and bound is in place but requires further development.

**Known Limitations**
- Branch and bound implementation is incomplete
- Limited error handling and input validation
- Basic output formatting
- No sensitivity analysis

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
