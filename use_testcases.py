import expected_distance as exp_dist
from result_analysis import analys_data, requirement_check, confidence_interval

cases = ['W', 'X', 'Y', 'Z']

str_models = ['WAG', 'LG', 'PAM', 'JTT']

seq_gen_dist = ["('A', 'B') is  0.80, ('A', 'C') is  1.0, ('B', 'C') is 1.2",
"('A', 'B') is  0.004, ('A', 'C') is  0.301, ('B', 'C') is  0.302",
"('A', 'B') is  0.301, ('A', 'C') is  0.701, ('B', 'C') is  1.0",
"('A', 'B') is  2.4, ('A', 'C') is  2.7, ('B', 'C') is 2.9"]



def all_testcases():
    """Finds the expected distance for all trees, all sets of sequences (1,2,3)
       with all models in the folder ./testcases. """
    for tree in range(len(cases)):
        run_testcase(tree)


def run_testcase(tree):
    """Takes a tree ('W', 'X', 'Y' or 'Z') as input and finds the
       expected distance for the tree for all sets of sequences (1,2 and 3)
       with all models from the folder ./testcases """

    total_mean = [0,0,0]
    data = []
    x1 = 1
    x2 = 36

    for nr in range(x1,x2): #1,2,3
    #    print("Sequences from ../"+cases[tree]+"/"+str(nr))
        for model in range(len(exp_dist.sub_models)):
            filename = "./testcases/"+cases[tree]+"/"+str(nr)+"/m"+str_models[model]+".fasta"
            print("\n Reading file.. :", filename)
            exps, stds, max_llhs, names, dists = (exp_dist.expected_dist(filename, exp_dist.sub_models[model], True))
            total_mean = [total_mean[i]+exps[i] for i in range(len(total_mean))]
            data.append([exps, stds, max_llhs, names, dists])

        #print("Newick", seq_gen_dist[tree])
    if plot_tree:
        exp_dist.draw_tree("./testcases/"+cases[tree]+".tree")
    if plot_analysis:
        long_estim_data = []
        for nr in range(x1,x2):
            for model in range(len(exp_dist.sub_models)):
                filename = "./testcases/"+cases[tree]+"/"+str(nr)+"/m"+str_models[model]+".fasta"
                exps, stds, max_llhs, names, dists = (exp_dist.expected_dist(filename, exp_dist.sub_models[model], False))
                long_estim_data.append([exps, stds, max_llhs, names, dists])
    #    analys_data(data, total_mean, tree, long_estim_data)
        requirement_check(data, tree, long_estim_data)
    #    confidence_interval(data, total_mean)
    #    confidence_interval(long_estim_data, total_mean)
exp_dist.plot_likelihood_function = False # Set as True to plot the likelihood function
exp_dist.plot_exp_dist = False # Set as True to plot exp
plot_tree = False  #Set as True to plot Newick Tree
plot_analysis = False  #Set as True to plot analysis compared to slow estimator
exp_dist.plot_posterior = False  #Set as True to plot posterior distribution


run_testcase(0) # cases[0] -> 'W'
#run_testcase(1) # cases[1] -> 'X'
#run_testcase(2) # cases[2] -> 'Y'
#run_testcase(3) # cases[3] -> 'Z'
#all_testcases()
