import math
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input radius",type = int)
radius = parser.parse_args()
print("area =",2*(math.pi)*radius,"circumference =",
      (math.pi)*radius**2, sep = " ")