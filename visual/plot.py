import sys
import json
import matplotlib.pyplot as plt
import math

with open ("out.json") as file :
  data = json.load(file)

potential = [p["V"] for p in data]
hole_conc = [math.log(p["p"]) for p in data]
elec_conc = [math.log(p["n"]) for p in data]

if str(sys.argv[1]) == "potential" :
  plt.plot(potential, lw = 3, label = "V")
  plt.title("Potential Curve")
  plt.xlabel("Grid points")
  plt.ylabel("Voltage V")
  
elif str(sys.argv[1]) == "concentration" :
  plt.plot(hole_conc, lw = 3, color = "r", label = "p")
  plt.plot(elec_conc, lw = 3, color = "purple", label = "n")
  plt.title("Concentration of holes and electrons")
  plt.xlabel("Grid points")
  plt.ylabel("log(concentration)")

plt.legend()
plt.grid(True)
plt.show()