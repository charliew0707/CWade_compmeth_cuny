# Comp Methods Homework 2: Integration


"""
HW2 – Exercise 5.3: Integration
Consider the integral

E(x) = int^x_0 e^{-t^2} dt  

a) Write a program to calculate E(x) for values of x from 0 to 3 in steps of 0.1. 
Choose for yourself what method you will use for performing the integral 
and a suitable number of slices.
b) When you are convinced your program is working, 
extend it further to make a graph of E(x) as a function of x.
Note that there is no known way to perform this particular integral 
analytically, so numerical approaches are the only way forward.

"""
import numpy as np
import argparse
import matplotlib.pyplot as plt

# function to define the integrand

def f(t):
    return np.exp(-t**2)

# function to define the integral, start with trapezoidal

def trap_int(a,b,N):

    if a==b:
        return 0.0
    
    dt = (b-a)/N
    t_vals = np.linspace(a,b,N+1)
    y = f(t_vals)
    s = 0.5 * y[0] + 0.5 * y[-1] + np.sum(y[1:-1])

    return float((dt * s))

# Function that integrates but uses Simpsons method:

def simp_int(a,b,N):

    if a==b:
        return 0.0
    
    if N % 2 != 0:
        raise ValueError("N must be an even integer for Simpson's rule.")
    
    dt = (b-a)/N
    t_vals = np.linspace(a,b,N+1)
    y = f(t_vals)
    s = y[0] + y[-1] + 4*np.sum(y[1:-1:2])+2*np.sum(y[2:-1:2])

    return float((dt/3)*s)



# function to calculate E(x), which add dynamic bound to b:

def E_numpy(xmax,dx,N_slices,int_function):

    num_points = int(round(xmax / dx)) + 1

    x_values = np.linspace(0,xmax,num_points) # creates np array from 0 to input, with 31 steps between, gives points 0.0 to 3.0 by 0.1
    if int_function == "simp":
        y_values = np.array([simp_int(0,x,N_slices) for x in x_values])
    else:
         y_values = np.array([trap_int(0,x,N_slices) for x in x_values]) 
    return x_values,y_values


# create simple function to plot

# function to plot with details - used gemini to label plot faster 
def plot_xy(x_data, y_data, xmax, dx, N, method):
    plt.figure(figsize=(8, 5))
    plt.plot(x_data, y_data, "b-o", label="$E(x) = \\int_0^x e^{-t^2} dt$")

    # Title showing the function
    plt.title("Numerical Evaluation of $E(x) = \\int_0^x e^{-t^2} dt$", fontsize=12)
    plt.xlabel("x")
    plt.ylabel("E(x)")
    plt.grid(True, linestyle="--", alpha=0.6)

    # Info blurb on the plot
    info_text = (
        f"Method: {method.upper()}\n"
        f"Slices (N): {N}\n"
        f"Step (dx): {dx}\n"
        f"Max x: {xmax}"
    )
    # Places the text box in upper-left inside axes coordinates
    plt.gca().text(
        0.05,
        0.92,
        info_text,
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8),
    )

    plt.legend(loc="center right")
    plt.tight_layout()
    plt.show()






# call function, use argparse that defaults to problem values, input steps, xmax value, and increment of x
def main():

    parser = argparse.ArgumentParser(description="Calculate and plot E(x) via numerical integration.")
    
    # 2. Add your 3 command-line arguments (Fill in the options!)
    parser.add_argument("-N", type=int , default=100 , help="Enter number of integration steps" )
    parser.add_argument("-xmax", type=float , default=3.0 , help="Enter upper bound of x for E(x)" )
    parser.add_argument("-dx", type=float , default=0.1 , help="Enter step size for x increment" )
    parser.add_argument("-m", choices=["trap", "simp"], default="trap", help="trap or simp for integration technique")
    
    # 3. Parse the arguments from the terminal
    args = parser.parse_args()
    x_data,y_data = E_numpy(args.xmax, args.dx, args.N,args.m)

    #Print x and E(x) data to check numbers before plotting
    print(f"{'x':<5} | {'E(x)':<10}")
    print("-" * 18)
    for x, val in zip(x_data, y_data):
        print(f"{x:.1f}   | {val:.6f}")
    
    # Call plot function
    plot_xy(x_data, y_data, args.xmax, args.dx, args.N, args.m)
    

if __name__ == "__main__":
    main()


# PRACTICE AND UNUSED FUNCTIONS


# trapezoidal integration with for loop not np
"""
def trap_int_loop(a,b,N):

    if a==b:
        return 0.0

    h = (b-a)/N
    s = 0.5*f(a)+0.5*f(b)

    for k in range(1,N):
        s += f(a+h*k)

    return s * h

"""
# old function for E(x):

"""
def E():
   
   # for ease of graphing, save data as separate lists
    x_val = []
    y_val = []

    for i in np.arange(0,3.1,0.1):
        x = round(i,1)
        I = float(trap_int(0,x,100))
        
        x_val.append(x)
        y_val.append(I)
    
    return x_val,y_val

"""