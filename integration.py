import matplotlib.pyplot as plt

def integrate_values(rule, xvalues, yvalues):
    return rule(yvalues, [xvalues[0], xvalues[-1]], len(xvalues))


def trapezoidal(list, interval, pnr):
    """ The function takes a function as intergrand, an interval and a
        parition number as input and returns a numerical approximation of the
        integral using Trapezoidal's rule. """
    a = interval[0]
    b = interval[-1]

# partition the interval [a,b] into pnr equal subintervals
    deltax = (b-a)/pnr
    sum = list[0] #first subinterval
    for i in list[1:]:
         #iteration for subintervals between a and b
         #deltax is added until the last interval is reached.
        a += deltax
        sum += 2*i
    sum += list[pnr-1] #last subinterval
    return (deltax/2)*sum

def simpson(y_values, x_values, pnr):
    """ The function takes a list of values and returns a numerical
        approximation of the
        integral using Simpson's rule."""
    a = x_values[0]
    b = x_values[-1]


    deltax = (b-a)/pnr
    sum = y_values[0] #first subx_values
    n = 0 #count iterations
    for y in y_values[1:-1]:
        n += 1
        if n%2 == 0: #if n is even
            sum += 2*y
        else: #else --> n is odd
            sum += 4*y
    sum += y_values[-1] #last subx_values
    return (deltax/3)*sum

def plot_integration(x_values, y_values, name, ylabel):
    plt.plot(x_values,y_values, label = "f")
    plt.xlabel('distances')
    plt.ylabel(ylabel)
    plt.title(name)
    plt.legend(loc='best')
    plt.show()
