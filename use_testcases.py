import expected_distance as exp_dist

cases = ['W', 'X', 'Y', 'Z']

str_models = ['WAG', 'LG', 'PAM', 'BLOSUM', 'JTT', 'CPREV', 'MTART', 'MTREV']

seq_gen_dist = ["('A', 'B') is  0.80, ('A', 'C') is  1.0, ('B', 'C') is 1.2",
"('A', 'B') is  0.004, ('A', 'C') is  0.301, ('B', 'C') is  0.302",
"('A', 'B') is  0.301, ('A', 'C') is  0.701, ('B', 'C') is  1.0",
"('A', 'B') is  2.4, ('A', 'C') is  2.7, ('B', 'C') is 2.9"]


def all_testcases():
    """Finds the expected distance for all trees, all sets of sequences (1,2,3)
       with all models in the folder ./testcases. """
    for tree in range(len(cases)):
        for nr in range(1,4):

            for model in range(len(exp_dist.sub_models)):
                filename = "./testcases/"+cases[tree]+"/"+str(nr)+"/m"+str_models[model]+".fasta"
                print("\n Reading file.. :", filename)
                exp_dist.expected_dist(filename, exp_dist.sub_models[model])

            print("Newick", seq_gen_dist[tree])

        if plot_tree:
            exp_dist.draw_tree("./testcases/"+cases[tree]+".tree")


def run_testcase(tree):
    """Takes a tree ('W', 'X', 'Y' or 'Z') as input and finds the
       expected distance for the tree for all sets of sequences (1,2 and 3)
       with all models from the folder ./testcases """

    for nr in range(1,4): #1,2,3
        print("Sequences from ../"+cases[tree]+"/"+str(nr))
        for model in range(len(exp_dist.sub_models)):
            filename = "./testcases/"+cases[tree]+"/"+str(nr)+"/m"+str_models[model]+".fasta"
            print("\n Reading file.. :", filename)
            exp_dist.expected_dist(filename, exp_dist.sub_models[model])
        print("Newick", seq_gen_dist[tree])
    if plot_tree:
        exp_dist.draw_tree("./testcases/"+cases[tree]+".tree")

exp_dist.plot_likelihood_function = False # Set as True to plot the likelihood function
exp_dist.plot_exp_dist = False # Set as True to plot exp
plot_tree = True  #Set as True to plot Newick Tree


run_testcase(0) # cases[0] -> 'W'
#run_testcase(1) # cases[0] -> 'X'
#run_testcase(2) # cases[0] -> 'Y'
#run_testcase(3) # cases[0] -> 'Z'
#all_testcases()
