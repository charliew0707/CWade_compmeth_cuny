# Computationa Physics Homework 1: Freefall 

"""

Write a program to calculate the time it takes for a ball to drop 
from a user speficified height to reach the ground. Use argparse. 
Also allow user to choose different values of gravity.
Add any other features that may be interesting.

"""

from math import sqrt
import argparse

# To improve: add other keywords for different planets maybe use a dict, 
# add try and except to catch errors, zero division, unphysical inputs
# add drag to calculation

# Function to calculate time for ball to fall from tower:

def time_to_fall(h, g):
    return round(sqrt((2*h)/g),4)

def main():
    parser = argparse.ArgumentParser(description="Freefall Calculator")
    parser.add_argument("height", type=float, help="Enter height of tower")
    parser.add_argument("--gravity", type=float, default=9.8, help="Enter gravity of planet")
    args = parser.parse_args()
    print(f"{time_to_fall(args.height, args.gravity)} seconds")

if __name__ == "__main__":
    main()