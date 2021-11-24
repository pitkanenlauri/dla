# Diffusion-limited aggregation

Diffusion-limited aggregation simulation written in Python.

Diffusion-limited aggregation (DLA) is the process whereby particles undergoing a random walk due to Brownian motion cluster together to form aggregates of such particles.

This program simulates such process with chosen grid size and particle count creating an image of the formed cluster. You can also calculate the fractal dimension of the created cluster.

## How to use

**Requirements:** Python 3.x.x, NumPy, SciPy and Matplotlib

**How to run:**
- If you want to **create cluster** of x times x grid with y particles, for example with values x=200 and y=2000: 
  - run with command: python dla.py 200 2000
  - -> cluster.txt file is created in the same directory with cluster xy-point pairs in separate rows.
- For **calculating fractal dimension of cluster**:
  - run with command: python dla.py name_of_the_cluster_file.txt
  - -> Shows cluster, the fit to the curve and prints out the estimated value for fractal dimension.

- **NOTICE!** As this is done in Python and also because of the nature of the simulation it can take a while to run the program with large number of particles. (many thousands or tens of thousands or more)

## Screenshots

