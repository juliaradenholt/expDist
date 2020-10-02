import matplotlib.pyplot as plt
import numpy as np

def integrate_values(rule, x_values, y_values, std = False):
    return rule(y_values, x_values, std)


def trapezoidal(y_values, x_values, bool_std):
    """ The function takes a function as intergrand, an interval and a
        parition number as input and returns a numerical approximation of the
        integral using Trapezoidal's rule. """
    a = x_values[0]
    b = x_values[-1]
    y_mul_dist = np.multiply(interval, len(x_values))
# partition the interval [a,b] into len(x_values) equal subintervals
    deltax = (b-a)/len(x_values)
    sum = y_values[0] #first subinterval
    std = y_values[0]**2
    for i in y_values[1:]:
         #iteration for subintervals between a and b
         #deltax is added until the last interval is reached.
        a += deltax
        sum += 2*i
        std += 2*(i**2)
    sum += y_values[len(x_values)-1] #last subinterval
    std += y_values[len(x_values)-1]**2
    if bool_std:
        return (deltax/2)*sum, (deltax/2)*std
    else:
        return (deltax/2)*sum

def simpson(y_values, x_values, bool_std):
    """ The function takes a y_values of values and returns a numerical
        approximation of the
        integral using Simpson's rule."""
    a = x_values[0]
    b = x_values[-1]
    y_mul_dist = np.multiply(x_values, y_values)

    deltax = (b-a)/len(x_values)
    sum = y_values[0] #first subx_values
    std = y_mul_dist[0]
    #print(sum, std)
    #std_lst = [y_values[0]**2]
    n = 0 #count iterations
    for y in y_values[1:-1]:
        n += 1
        if n%2 == 0: #if n is even
            sum += 2*y
            std += 2*y_mul_dist[n]
        #    std_lst.append(2*(y**2))
        else: #else --> n is odd
            sum += 4*y
            std += 4*y_mul_dist[n]
        #    std_lst.append(4*(y**2))
        #print(std)
    sum += y_values[-1] #last subx_values
    std += y_mul_dist[-1]
    if bool_std:
        return (deltax/3)*sum, (deltax/3)*std
    else:
        return (deltax/3)*sum

def plot_integration(x_values, y_values, name, ylabel):
    plt.plot(x_values,y_values, label = "f")
    plt.xlabel('distances')
    plt.ylabel(ylabel)
    plt.title(name)
    plt.legend(loc='best')
    plt.show()
