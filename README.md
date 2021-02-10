# ExpDist


This is my bachelor's thesis project.
The project was to implement an efficient estimator for the expected value of evolutionary dsitances.

Title: ExpDist, an estimator for expected values of evolutionary distances
University: Stockholm University

## Setup 
Clone the current version of the project.

```
git clone https://github.com/juliaradenholt/expDist.git
```

Before the program can be configured, we must install some libraries without these the program will not run

```
$ pip install numpy
$ pip install biopython
```
Optional library, which has to be installed if the user wants to generate plots. 
```
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

## run testcases
In the folder called test_cases, four phylogentic trees are available. Sequences were generated from these trees using PAM, WAG, LG and JTT as substitutions models. The sequences were generated 35 times per tree and substitution model. Therefore a total of 420 alignments per model are available. The substitutions model which was used to generate the sequences is also used as input for ExpDist.
To test the estimator with testcase 'W', run

```
python3 use_testcases.py
```
To test the estimator with the rest of the test cases, uncomment the following parts in use_testcases.py

```
#run_testcase(1) # cases[1] -> 'X'
#run_testcase(2) # cases[2] -> 'Y'
#run_testcase(3) # cases[3] -> 'Z'
#all_testcases() 
```
To see plots of the tree, the expected value, the posterior distribution or the likelihood,
set the following variables in use_testcases.py as true:
```
exp_dist.plot_likelihood_function = False # Set as True to plot the likelihood function
exp_dist.plot_exp_dist = False # Set as True to plot exp
plot_tree = False  #Set as True to plot Newick Tree
exp_dist.plot_posterior = False  #Set as True to plot posterior distribution
```
To receieve analysis of the estimator, set the following variable in use_testcases.py as true:
```
plot_analysis = False  #Set as True to plot analysis compared to slow estimator
```


