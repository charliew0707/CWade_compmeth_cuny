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

# Adding Planets
PLANETS = {
    "mercury": 3.7,
    "venus"  : 8.87,
    "earth"  : 9.81,
    "moon"   : 1.62,
    "mars"   : 3.72,
    "jupiter": 24.79,
    "saturn" : 10.44,
    "uranus" : 8.69,
    "neptune": 11.15,
}

# Function to calculate time for ball to fall from tower:

def time_to_fall(h, g):
    return round(sqrt((2*h)/g),4)

def main():
    
    parser = argparse.ArgumentParser(description="Freefall Calculator")
    parser.add_argument("height", type=float, help="Height of tower in meters")
    parser.add_argument("--gravity", type=float, default=None,
                        help="Custom gravity value (m/s^2)")
    parser.add_argument("--planet", type=str, default="earth",
                        choices=PLANETS.keys(),
                        help="Planet name (default: earth)")
    args = parser.parse_args()

    # --gravity overrides --planet if both are given
    g = args.gravity if args.gravity is not None else PLANETS[args.planet]

    print(f"Height  : {args.height} m")
    print(f"Gravity : {g} m/s^2")
    print(f"Time    : {time_to_fall(args.height, g)} seconds")

if __name__ == "__main__":
    main()