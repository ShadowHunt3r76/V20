# Linear Programming and Integer Programming Solver

A .NET 9 application that solves Linear Programming (LP) and Integer Programming (IP) models, including the Knapsack problem. This project was developed as part of the LPR 381 (Linear Programming) course, focusing on implementing core optimization algorithms with an emphasis on educational value and code clarity.

## 🚀 Getting Started

### Prerequisites
- [.NET 9.0 SDK](https://dotnet.microsoft.com/download/dotnet/9.0) or later
- (Optional) Visual Studio 2022 or VS Code with C# extension

### Installation
1. Clone the repository
2. Navigate to the project directory:
   ```bash
   cd V20
   ```

### Running the Application

#### Method 1: Using Command Line Arguments
```bash
dotnet run --project LinearProgrammingApp -- [input_file] [output_file]
```

Example:
```bash
dotnet run --project LinearProgrammingApp -- knapsack_problem.txt solution.txt
```

#### Method 2: Interactive Mode
1. Run without arguments:
   ```bash
   dotnet run --project LinearProgrammingApp
   ```
2. Follow the on-screen prompts to:
   - Enter the path to your input file
   - Select a solver from the menu
   - View and save the solution

## 🧩 Supported Problems

### 1. Linear Programming (LP)
- **Primal Simplex Algorithm**
- **Revised Primal Simplex Algorithm**

### 2. Integer Programming (IP)
- **Branch and Bound**
- **Cutting Plane Method**

### 3. Knapsack Problem
- **Dynamic Programming Solution**
- **Binary Decision Variables**

## 📝 Input File Format

The input file must follow this structure:

```
[objective] [coefficient1] [coefficient2] ... [coefficientN]
[constraint1_coefficients] [relation] [rhs]
[constraint2_coefficients] [relation] [rhs]
...
[var1_type] [var2_type] ... [varN_type]
```

### Example: Knapsack Problem
```
max +2 +3 +3 +5 +2 +4
+11 +8 +6 +14 +10 +10 <= 40
bin bin bin bin bin bin
```

### Example: Linear Program
```
max +3 +2
+1 +2 <= 6
+1 +1 <= 4
+ +1 <= 3
+ +
```

## 🛠️ Solvers

1. **Primal Simplex (LP)** - Best for small to medium linear programs
2. **Revised Primal Simplex (LP)** - More efficient for larger problems
3. **Branch and Bound (IP)** - For integer programming problems
4. **Cutting Plane (IP)** - Alternative method for IP problems
5. **Knapsack Solver** - Specialized for binary knapsack problems

## 📂 Project Structure

- `LinearProgrammingApp/` - Main console application
- `LinearProgramming.Parsing/` - Input parsing and model representation
- `LinearProgramming.Algorithms/` - Core optimization algorithms
  - `PrimalSimplex/` - Tableau-based simplex implementation
  - `BranchAndBound/` - Integer programming solver
  - `CuttingPlane/` - Alternative IP solver
  - `Knapsack/` - Specialized knapsack solver
- `LinearProgramming.Tests/` - Unit tests

## 📊 Output

Solutions are saved to the specified output file and include:
- Optimal solution (if found)
- Objective value
- Variable values
- Solution status (Optimal/Infeasible/Unbounded)

## 🧪 Testing

Run all tests with:
```bash
dotnet test
```

## 📚 Documentation

- XML documentation is available in the code
- Algorithm implementations include detailed comments
- Test cases serve as usage examples

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

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
