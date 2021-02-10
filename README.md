# ExpDist


This is my bachelor's thesis project.
The project was to implement an efficient estimator for the expected value of evolutionary dsitances.

Title: An efficient estimator for expected values of evolutionary distances.
University: Stockholm University

## Setup 
Clone the current version of the project.

```
git clone https://github.com/juliaradenholt/expDist.git
```

Before the program can be configured, we must install some libraries without which the program will not run

```
$ pip install numpy
$ pip install biopython
$ pip install mathlibplot
```
## Estimate expected values 

In the terminal, run 
```
python3 run_expdist.py
```

or in the Python Shell run:
```
import run_expdist.py
```

This will allow you to enter a substiution model and a filepath.
The program will then estimate and print the expected values of the evolutionary distances between the sequences in the file.
