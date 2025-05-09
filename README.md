# Assignment 2 - Numerical Linear Algebra

This repository contains solutions to the second assignment of the [Numerical Linear Algebra course at FGV/EMAp](https://emap.fgv.br/en/discipline/numerical-linear-algebra).

## Assignment Description

In this project, we tackle the least squares linear and polynomial regression problems by:
- Solving the systems using normal equations, as well as QR and SVD decompositions.
- Analyzing the condition number of the matrices involved to evaluate numerical stability.
- Discussing the sensitivity of the solutions in response to data perturbations.
- Generating plots to visualize and compare the outcomes of the implemented methods.

The problems are based on 100 equally spaced points in the interval [0, 1], and another dataset centered around the origin. We analyze the difference in conditioning.

## Features

- Implementation of least squares regression using multiple numerical methods.
- Detailed analysis of condition numbers and sensitivity to matrix perturbations.
- Generation of plots for comparative visualization.
- Modular and documented codebase for easy understanding and extensibility.

## Project Structure

- **assignment.ipynb**: Contains the source code for regression methods and numerical analyses.
- **assignment.typ**: Contains the Typst source code for the report
- **assignment.pdf**: The compiled report in PDF format.
- **plots/**: Folder with the plots generated in assignment.ipynb.
- **README.md**: This file, which explains the project.
- **LICENSE**: Contains the licensing details of this repository.


### Requirements

We recommend using a virtual environment to manage dependencies. You can create one using `venv`:

If you are using a Linux-based OS:

```bash
python -m venv coolvenv
source coolvenv/bin/activate
```

If you are using Windows:

- Immediately uninstall Windows and install _any_ Linux-based OS. 

Then, install the required packages using `pip`:

```bash

pip install -r requirements.txt
```

We use a Jupyter notebook for the implementation, which requires the following Python libraries:

- [`numpy`](https://numpy.org/): For numerical computations and matrix operations.
- [`matplotlib`](https://matplotlib.org/): For plotting and visualizing data.

You can install the required packages with:

```bash
pip install -r requirements.txt
```

## How to Run

After installing the required packages, you can run all cells in the Jupyter notebook `assignment.ipynb` to execute the regression analyses and generate the plots.

### Instructions

1. **Clone the repository:**

     ```bash
     git clone https://github.com/arthurabello/nla-assignment-2.git
     ```

2. **Navigate to the project directory:**

     ```bash
     cd nla-assignment-2
     ```

3. **Run the Jupyter Notebook (All Cells)**

4. **Read the report:**

   - The report is written in Typst. Read the PDF version `assignment.pdf` for a detailed explanation of the methods and results.

## License

The project is licensed under the terms specified in the [LICENSE](https://github.com/arthurabello/nla-assignment-2/blob/main/LICENSE) file.
