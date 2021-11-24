# /usr/bin/env python
import sys
import copy
from math import *
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def count_particles(center, boxsize, grid):
    """
    Count the number of points in the grid
    with the value of 1 (or True) for which
    the x and y indices differ from the coordinates
    [x0 ,y0] given in 'center' by at most boxsize.
    
    That is, count how many of the grid points in
    [x0-boxsize, x0+boxsize] x [y0-boxsize, y0+boxsize]
    are 1.
    """
    minx = center[0]-boxsize
    maxx = center[0]+boxsize
    miny = center[1]-boxsize
    maxy = center[1]+boxsize

    particles = 0
    for i in range(minx, maxx):
        for j in range(miny, maxx):
            if grid[i,j]:
                particles += 1
    return particles
    

def random_step(x,y):
    """
    Take a random step of length 1
    from the given coordinates (x,y)
    and return the new coordinates.
    """
    rnd = random.rand()
    if rnd < 0.25:
        return [x+1, y]
    elif rnd < 0.5:
        return [x-1, y]
    elif rnd < 0.75:
        return [x, y+1]
    else:
        return [x, y-1]
    return [x, y]


def neighbor_is_in_cluster(x,y, grid):
    """
    Check if any of the neighbors of the given point
    (x, y) are already a part of the cluster.
    A point is part of the cluster if that point
    has the value of 1 in the given array 'grid'.
    """
    try:
        if grid[x-1,y]:
            return True
        elif grid[x+1,y]:
            return True
        elif grid[x,y-1]:
            return True
        elif grid[x,y+1]:
            return True
        else:
            return False
    except:
        return False


def pick_random_point(center, distance):
    """
    Choose a random point at the given distance from
    the given center. Distance should be a real number
    and center should be an array of two real numbers,
    [x, y].
    """
    angle = random.rand()*2.0*pi    
    point = center + distance * np.array([cos(angle), sin(angle)])
    return np.array( [int(point[0]), int(point[1])] )


def distance_from_center(point, center):
    """
    Returns the integer value of distance of a given point from center.
    Point and center should be arrays of two real numbers, [x, y].
    """
    x_1, y_1 = center[0], center[1]
    x_2, y_2 = point[0], point[1]

    # Use Pythagorean theorem to calculate the distance. :)
    return ceil(sqrt( (x_2-x_1)**2 + (y_2-y_1)**2) )


def diffusion_simulation(size, particles):
    """
    Runs a diffusion simulation.

    Parameters:
    # size, determines a 'size' x 'size' grid from which center
    point the cluster will be constructed.
    # particles, sets the amount of particles the final cluster
    will consist.

    Function:
    Particles are placed some distance from cluster radius and
    are set to random walk until they reach the cluster. When
    reached they attach to it.

    Returns:
    # grid, matrix consisting cluster as its location defined by
    elements with value 1 with the rest being 0.
    # center, the center of the grid as an array, (x, y).
    """
    # Clusters longest arm distance from center i.e. cluster radius.
    r_cluster = 1

    # Limit walkers wandering too far from cluster.
    max_distance = 5
    
    # Create a grid of integers of the size 'size' x 'size'
    # with all points having value 0.
    grid = np.zeros((size, size), dtype=int)

    # Set the center of the grid.
    center = [size//2, size//2]
    
    # Assign center point of the grid to have value 1 as a starting point.
    grid[center[0]][center[1]] = 1

    # Add the parameter 'particles' times particles to cluster.
    for i in range(particles):

        # Put particle near cluster.
        location = pick_random_point(center, r_cluster)
        
        # Do the random walk.
        walk = True
        while walk:
            # Change location by taking the random step.
            location = random_step(location[0], location[1])

            # Keep track of how far the particle is from the center.
            d = distance_from_center(location, center)

            # See if particle finds cluster.
            if neighbor_is_in_cluster(location[0], location[1], grid):
                walk = False

                # Add then particle to the cluster.
                grid[location[0]][location[1]] = 1
                
                # Update r_cluster if new the point is farthest one so far.
                if d > r_cluster:
                    r_cluster = d
            
            # Keep particles from reaching too far from r_cluster.
            if d > r_cluster + max_distance:
                location = pick_random_point(center, r_cluster)

    return grid, center


def read_grid_from_file(filename):
    """
    Read a cluster from a file.
    
    The first line of the file should contain the
    coordinates of the initial seed of the cluster.
    This should be followed by the coordinates of all
    other points in the cluster, line by line:
    
    x0 y0
    x1 y1
    x2 y2
    ...
    
    The function returns a grid describing the cluster and the
    position of the initial seed.
    """
    f = open(filename)
    lines = f.readlines()
    f.close()
    maxindex = 0
    
    parts = lines[0].split()
    center = [ int(parts[0]), int(parts[1]) ]
        
    grid = np.array( [ [False]*(center[0]*2) ]*(center[1]*2) )    

    for line in lines:
        parts = line.split()
        if len(parts) > 0:
            grid[ int(parts[0]), int(parts[1]) ] = True

    return grid, center


def linearfit(x, a, b):
    """
    The function f(x) = ax + b
    """
    return a * x + b


def calculate_fractal_dimension(grid, center):
    """
    Calculate the fractal dimension of the cluster.
    
    The input:
    * grid: an array of points where the value of 1 denotes 
            a point that belongs in the cluster
    * center: the coordinates of the initial seed of the cluster    
    """
    # store the statistics in these lists
    particlecount = []
    logsize = []
    logcount = []
    R = 1 # box size

    # Loop over increasing box area.
    walk = True
    while walk:
        # Calculate how many particles are found for box sized R.
        N = count_particles(center, R, grid)

        # Break out of the loop if increasing R doesn't increase N,
        # meaning we've reached the outermost edge of the cluster.
        if N > 10 and N == particlecount[-1]:
            break
        
        particlecount.append(N)

        # Calculate logarithms ln N and ln R and store values.
        logsize.append(log(R))
        logcount.append(log(N))

        R += 1
        
    # fit a straight line to the data
    logsize = np.array(logsize)
    logcount = np.array(logcount)
    
    # The first and last values in the data to be included in the fit.
    # Eyeballed to represend the linear part of the fit. ;)
    first_value = 0
    last_value = -160

    popt, pcov = curve_fit(linearfit, 
                            logsize[first_value:last_value],
                            logcount[first_value:last_value], 
                            p0 = [1.5, 2.0]) # initial guess for parameters
    print( "Estimated value for the fractal dimension: " + str(popt[0]) )


    print(logsize)
    plt.plot(logsize, logcount)
    plt.plot(logsize, linearfit(logsize, popt[0], popt[1]))
    plt.show()


def show_cluster(grid):
    """
    Plot the cluster.
    """
    plt.imshow(grid, cmap='Greys', interpolation='nearest')
    plt.show()


def write_cluster_file(grid, center):
    """
    Write the cluster in a datafile which can be
    later read in using read_grid_from_file().
    """
    size = len(grid[0])
    writelines = str(center[0])+" "+str(center[0])+"\n"
    
    for i in range(size):
        for j in range(size):
            if grid[i,j]:
                writelines += str(i)+" "+str(j)+"\n "
    f = open('cluster.txt','w')
    f.write(writelines)
    f.close()
    
    

def main(args):
    """
    The main program.
    
    If you wish to read in the results of a previous
    simulation, give the name of the file to read in
    'filename'
    """
    grid = []
    center = []
    filename = args[0] # name of the file to read
    try:
        # read a file containing cluster data
        grid, center = read_grid_from_file(filename)
    except:
        # if no valid file was found, run the simulation
        size = int(args[0])
        particles = int(args[1])
        grid, center = diffusion_simulation(size, particles)
        
    show_cluster(grid)
    
    calculate_fractal_dimension(grid, center)

    write_cluster_file(grid, center)


#
# Run the main program if this file is executed
# but not if it is read in as a module.
#
if __name__ == "__main__":
    random.seed()
    main(sys.argv[1:])
