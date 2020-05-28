import matplotlib.pyplot as plt
import numpy as np


true_values = [[0.80, 1.0, 1.2],[0.004,0.301,0.302],[0.302,0.701,1.0],[2.4,2.7,2.9]]


def confidence_interval(data, total_mean):
    #print(len(data[0][0]))
    #print(total_mean)
    #print(len(data))

    plt.figure()
    for d in data:
        exps, stds, max_llhs, names, dists = d
        n = len(data)
        mean = [i/n for i in total_mean]
        margin_error = []

        calculation_num = [str(i) for i in names]
        plt.plot(exps, calculation_num, '.', color = 'grey', markersize = 2)


    plt.title("Sample of expected values")
    plt.show()





def requirement_check(data, tree, long_estim_data):
    total = 0
    n = 0
    for i in range(len(data)):
        fast_exps, fast_stds, fast_max_llhs, names, fast_dists = data[i]
        long_exps, long_stds, long_max_llhs, names, long_dists = long_estim_data[i]

        for j in range(len(names)):
            total +=1
            if ((fast_exps[j]/long_exps[j] < 1.1) and (fast_exps[j]/long_exps[j] > 0.9)) or ((true_values[tree][j] < (fast_exps[j]+fast_stds[j])) and (true_values[tree][j] > (fast_exps[j]-fast_stds[j]))):
                n+=1
            else:
                print("No requirement reached..")
                print("Long est:", long_exps[j])
                print("Fast est:", fast_exps[j])
                print(fast_exps[j]/long_exps[j])
                print(fast_exps[j]-fast_stds[j], true_values[tree][j], fast_exps[j]+fast_stds[j])

    print("Requirement reached", str(n)+"/"+str(total), "times")

def analys_data(data, total_mean, tree, long_estim_data = []):
    n = 0
    total = 0
    if long_estim_data == []:
        for i in data:

            exps, stds, max_llhs, names, dists = i

            interval = [[i[1][0],i[1][-1]]  for i in dists]
            calculation_num = [str(i) for i in names]
            for i in range(len(calculation_num)):
                    plt.plot(interval[i], [calculation_num[i],calculation_num[i]], ':', color = 'grey')
            plt.plot([], [], ':', color = 'grey', label = 'Interval (Discretization)')
            plt.plot(max_llhs, calculation_num, 'ro', color = 'red', label = "Maximum likelihood")
            plt.plot(true_values[tree], calculation_num, 'ro', color = 'green', label = "Newick Tree values")
            plt.errorbar(exps, calculation_num, None, stds, color = 'black', label="Expected value with standard deviation", linestyle='None', markersize=5, marker='^')
            plt.legend(loc='best')
            plt.show()
    else:
        for i in range(len(data)):

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12,6))
            fig.suptitle('Analysis')
            ax1.title.set_text('Fast estimator (ExpDist)')
            ax2.title.set_text('Slow estimator')


            fast_exps, fast_stds, fast_max_llhs, names, fast_dists = data[i]
            long_exps, long_stds, long_max_llhs, names, long_dists = long_estim_data[i]

            fast_interval = [[i[1][0],i[1][-1]]  for i in fast_dists]
            long_interval = [[i[1][0],i[1][-1]]  for i in long_dists]
            calculation_num = [str(i) for i in names]

            #print("Discretization points", len(fast_dists[0][1]),len(fast_dists[1][1]),len(fast_dists[2][1]))

            for i in range(len(calculation_num)):
                    ax1.plot(long_exps[i], calculation_num[i], 'v', color = 'orange', markersize = 10)
                    ax1.plot(fast_interval[i], [calculation_num[i],calculation_num[i]], ':', color = 'grey')
                    ax2.plot([0,3.0], [calculation_num[i],calculation_num[i]], ':', color = 'grey')


            ax1.plot([], [], ':', color = 'grey', label = 'Interval (Discretization)')
            ax2.plot([], [], ':', color = 'grey', label = 'Interval (Discretization)')
            ax1.plot(fast_max_llhs, calculation_num, 'ro', color = 'red', label = "Maximum likelihood")

            ax1.plot(true_values[tree], calculation_num, 'ro', color = 'green', label = "Newick Tree distance")

            ax1.errorbar(fast_exps, calculation_num, None, fast_stds, color = 'black', label="Expected value with\n standard deviation", linestyle='None', markersize=5, marker='^')

            ax2.plot(long_max_llhs, calculation_num, 'ro', color = 'red', label = "Maximum likelihood")
            ax2.plot(true_values[tree], calculation_num, 'ro', color = 'green', label = "Newick Tree distance")
            ax2.errorbar(long_exps, calculation_num, None, fast_stds, color = 'black', label="Expected value with\n standard deviation", linestyle='None', markersize=5, marker='^')
            ax1.plot([], [],'v', color = 'orange', label = "Mean (Slow estimator)")

            info_string = "Pairwise distance between "+str(names)+"\n"
            newick_string = "Newick Tree distance: "+str(true_values[tree])+"\n\n"
            fast_data_string = "Expected value for distance: "+str((round(fast_exps[0],2), round(fast_exps[1],2), round(fast_exps[2],2)))+"\nMaximum likelihood: "+str( (round(fast_max_llhs[0],2), round(fast_max_llhs[1],2), round(fast_max_llhs[2],2)))
            long_data_string = "Expected value for distance: "+str((round(long_exps[0],2), round(long_exps[1],2), round(long_exps[2],2)))+"\nMaximum likelihood: "+str((round(long_max_llhs[0],2), round(long_max_llhs[1],2), round(long_max_llhs[2],2)))
            difference_string = "\n\nPercentage difference between \nExpDist and Slow estimator: "+str( (round(np.abs((fast_exps[0]/long_exps[0])-1),3),round(np.abs((fast_exps[1]/long_exps[1])-1),3),round(np.abs((fast_exps[2]/long_exps[2])-1),3)))
            ax1.text(0.02, 0.5, info_string+newick_string+"ExpDist: \n"+fast_data_string+"\n\nSlow estimation: \n"+long_data_string+difference_string, fontsize=12, transform=plt.gcf().transFigure)

            plt.subplots_adjust(left=0.5)
            ax1.legend(loc='upper left', bbox_to_anchor=(0.8,0.8), fontsize='small')
            ax2.legend(loc='upper left', bbox_to_anchor=(0.8,0.8), fontsize='small')
            plt.show()
    #print(np.multiply(total_mean,1/len(data)))
