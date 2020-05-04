import numpy as np
from initialize_data import readFasta, initializeSequenceData, draw_tree, sub_models
from integration import trapezoidal, simpson, plot_integration, integrate_values
from text import color

integration_method = simpson
plot_likelihood_function = False
plot_exp_dist = False
print_maximum_likelihood = False



def replacement_matrix(model, t):
    """ Returns P(t) = e^(Qt) """
    return model.get_replacement_probs(t)



def prior_prob(model, distances):
    """ Takes a model and a set of distances and
        returns a list of replacement matrices. """
    repl_matrices = []
    for d in distances:
        repl_matrices.append(replacement_matrix(model, d))
    return repl_matrices



def likelihood_estimation(alignmentMatrix, repl_matrices, equilibrium_freq):
    """ Takes the difference position matrix for two sequences and
        a list of replacement matrices. Returns a list of values
        correspoding to the likelihood. """
    likelihood = []
    freqs = (np.prod(equilibrium_freq)) #product of frequencies πa1 * πa2 * ...

    for m in repl_matrices:
        mul = np.power(np.array(m), alignmentMatrix)
        lh = np.prod(mul)*freqs # Pr(a1)^x *  πa1 * Pr(a2)^y * πa2 * ...
        likelihood.append(lh)

    return likelihood




def expected_dist(filename, model):
    read_align = readFasta(filename, model)
    all_alignments = [initializeSequenceData(read_align[0][i[0]], read_align[0][i[1]], model) for i in read_align[2]]
    names = read_align[1]
    length = len(read_align[0][0])
    print("\nEstimating sequences with length: "+color.BOLD+str(length)+color.END+" with model: "+color.BOLD+str(model.get_name())+color.END+"\n")

    for i in range(len(all_alignments)):

        aa_map = all_alignments[i][0]

        #distances = [x for x in np.arange(0, 3, 0.005)]
        distances = all_alignments[i][1]
        #print("Number of discretization points", len(distances))

        freq_appearence = all_alignments[i][2] # [π_a,  π_r, ...]

        repl_matrices = prior_prob(model, distances) # [e^Qd1, e^Qd2, ...]
        likelihood = likelihood_estimation(aa_map, repl_matrices, freq_appearence) # [Pr(a,b|d1), Pr(a,b|d2), ..]
        total_llh = (integrate_values(integration_method, distances, likelihood)) # ∑ Pr(a,b|d)

        posterior_dist = np.multiply(distances, likelihood) # [d1*Pr(d1|a,b), d2*Pr(d2|a,b), ...]

        expected_distance = (1/total_llh)*(integrate_values(integration_method, distances, posterior_dist))
        print("Expected distance for "+color.BOLD+str(names[i])+color.END+" is ",
        color.BOLD+str(expected_distance)+color.END)


        if plot_likelihood_function:
            plot_integration(distances, np.multiply(likelihood, 1/total_llh), "Distances: "+str(names[i]), 'Pr(d|a,b)')

        if print_maximum_likelihood:
            print("ML/MAP = argmax(P(d|a,b)) =", distances[likelihood.index(np.max(likelihood))])

        if plot_exp_dist:
            plot_integration(distances, np.multiply(posterior_dist, 1/total_llh), "test", 'd*Pr(d|a,b)')
