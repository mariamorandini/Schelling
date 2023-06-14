import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import seaborn as sns
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy.signal import convolve2d
import warnings
from multiprocessing import Pool, cpu_count
import numpy as np
from scipy.ndimage import label

warnings.filterwarnings("ignore")


def step(n_steps, occupied_sites):
    
    """ This function returns the actual number of required steps to run the simulations, by definition:
    one step = one move per occupied site on average 
    occupied sites -- self.population
    """
    return n_steps*occupied_sites 


"""
class UnionFind:
    def __init__(self, max_labels):
        self.labels = [0] * max_labels
        self.labels[0] = 0
        self.n_labels = max_labels
    def find(self, x):
        y = x
        while self.labels[y] != y:
            y = self.labels[y]
        while self.labels[x] != x:
            z = self.labels[x]
            self.labels[x] = y
            x = z
        return y
    def union(self, x, y):
        self.labels[self.find(x)] = self.find(y)
        return self.find(x)
    def make_set(self):
        self.labels[0] += 1
        assert self.labels[0] < self.n_labels
        self.labels[self.labels[0]] = self.labels[0]
        return self.labels[0]


def hoshen_kopelman(matrix, m, n):
      uf = UnionFind(m * n // 2)
      for i in range(m):
          for j in range(n):
              if matrix[i][j]:
                  up = matrix[i - 1][j] if i > 0 else 0
                  left = matrix[i][j - 1] if j > 0 else 0
                  if up == 0 and left == 0:
                      matrix[i][j] = uf.make_set()
                  elif up == 0 or left == 0:
                      matrix[i][j] = max(up, left)
                  else:
                      matrix[i][j] = uf.union(up, left)
      return matrix"""

def dfs(row, col, cluster_type):
    """
    Helper function to perform depth-first search to find the size of the
    nearest neighbor cluster of the given type.
    """
    if row < 0 or row == n_rows or col < 0 or col == n_cols:
        return 0

    if visited[row, col] or city[row, col] != cluster_type:
        return 0

    visited[row, col] = True
    size = 1

    size += dfs(row - 1, col, cluster_type)
    size += dfs(row + 1, col, cluster_type)
    size += dfs(row, col - 1, cluster_type)
    size += dfs(row, col + 1, cluster_type)
    return size
    





##############
class Schelling:

    def __init__(self, size, tolerance, p, p_plus=.5, matrix=None):
        """This class is used to store all the datas of our model.
        If not previously specified it initializes a matrix representing our system,
        with the features given.
        It allows to greate an instance of tis class with a given initial condition
        if the matrix is different from m.
        """
        self.matrix = matrix
        self.T = tolerance
        self.p = p  # vacancy density
        self.p_plus = (1 - self.p) * p_plus
        self.p_minus = 1 - self.p_plus - self.p
        self.size = size

        # this allow to introduce the same initial condition to a system,
        # if not specified previously we initialize a matrix according to the given parameters
        #### TO BE ADDED --- checks that our model initial parameters match

        if np.array(self.matrix).any() != None:
            self.city = np.array(matrix)
        else:
            self.city = np.array(
                np.random.choice([1, -1, 0], size=self.size ** 2, p=[self.p_plus, self.p_minus, self.p]).reshape(
                    self.size, self.size))

        # population exact values
        self.vacancies = (self.city == 0).sum()  # Number of vacancies
        self.plus = (self.city == 1).sum()  # Number of agents of cluster 1
        self.minus = (self.city == -1).sum()  # Number of agents of cluster -1
        self.population = self.plus + self.minus  # Number of agents

        # internal variables
        self.kernel = np.ones((3, 3))  # used to compute Moore  neighborhood
        self.kernel[1, 1] = 0
        self.K = 2 * self.T - 1

    def get_location_people(self):
        """ Return a np.array of tuples of coordinates of people's positions """

        people = []

        for x in np.arange(self.size):
            for y in np.arange(self.size):
                if self.city[x, y] != 0:
                    people.append((x, y))

        return people

    def dissatisfied_number(self):
        """ Return the number of agents dissatisfied with their neighboor"""

        return len(self.get_dissatisfied_location())

    def satisfied_number(self):
        """ Return the number of agents satisfied with their neighboor"""
        return self.population - self.dissatisfied_number()
    


    def opposite_edges_ratio(self, matrix = True):
        """ returns the ratio of edges of opposite type, vs the number of edges between any agents """
        matrix = self.city
        kernel =  self.kernel
        agent_mask = np.abs(matrix)  # Mask for agent cells (1 or -1)
        agent_edges = convolve2d(agent_mask, kernel, mode='same', boundary='fill')  # Count agent edges

        opposite_mask = matrix < 0  # Mask for cells with opposite sign
        opposite_edges = convolve2d(opposite_mask, kernel, mode='same', boundary='fill')  # Count opposite edges

        opposite_edges_count = np.sum(opposite_edges)
        total_edges_count = np.sum(agent_edges)

        if total_edges_count == 0:
            return 0  # Avoid division by zero

        return opposite_edges_count / total_edges_count

    

    def get_dissatisfied_location(self, boundary="fill", fillvalue=0):
        """ Return the location of agents dissatisfied with their neighboor,
        if not specified, the boundary conditions imposed are periodic
        boundary = "wrap" """

        kws = dict(mode='same', boundary=boundary, fillvalue=fillvalue)  # periodicity of boundaries is included

        con_plus = convolve2d(self.city == 1, self.kernel, **kws)
        con_minus = convolve2d(self.city == -1, self.kernel, **kws)
        neighs = convolve2d(self.city != 0, self.kernel, **kws)

        plus_diss = (con_plus / neighs <= 1 - self.T) & (self.city == 1)
        minus_diss = (con_minus / neighs <= 1 - self.T) & (self.city == -1)

        dissatisfied = []

        for i in range(self.size):
            for j in range(self.size):
                if plus_diss[i, j]:
                    dissatisfied.append((i, j))
                if minus_diss[i, j]:
                    dissatisfied.append((i, j))


        return dissatisfied

    def run_dissatisfied(self, boundary="fill", fillvalue=0):
        
        """One time step of the simulation. Among the dissatisfied agents, one at random is chosen
        and moved to a satisfying location (i. e. where tolerance level is satisfied, not necessarly strictly).
        If no satisfying relocation or no dissatisfied agent is present, no change is done.
       """
        x = self.get_dissatisfied_location() 
        if len(x) == 0:
            # print("reached frozen configuration: No dissatisfied people are present")
            return None

        person_index = x[
            np.random.choice(np.arange(len(x)))]
        genre = self.city[person_index]
        

        # get the values ine each cell of the number of neighbors, and number of similar ones
        kws = dict(mode='same', boundary=boundary)  # periodicity of boundaries is included
        con = convolve2d((self.city == genre).astype(int), self.kernel, **kws)
        neighs = convolve2d((self.city != 0).astype(int), self.kernel, **kws)

        # new possible relocations are the ones that satisfy T condition and empty
        possible_relocations = ((con / neighs >= 1 - self.T)) & ((self.city == 0))
        
        # collects the indexes of the possible relocations
        
        index_rc = []  

        for x in np.arange(self.size):
            for y in np.arange(self.size):
                if possible_relocations[x, y]:
                    index_rc.append((x, y))

        if len(index_rc) == 0:
            # print("reached frozen configuration: no suitable relocation present")
            return None
        

        replace = index_rc[np.random.choice(np.arange(len(index_rc)))]
        self.city[person_index[0], person_index[1]] = 0
        self.city[replace] = genre

    def run_everybody(self, boundary="fill", fillvalue=0):
        
        """One time step of the simulation. Among ALL agents, one at random is chosen
        and moved to a satisfying location (i. e. where tolerance level is satisfied, not necessarly strictly),
        thus the utility does not necessarly increases at each time step.
        --- this adds randomness in our model, will be analyzed deeper later.
        If no satisfying relocation no change is done.
        The condition on the boundary is by default periodic. """
        # picking a random index 
        x = self.get_location_people()
        person_index = x[np.random.choice(np.arange(len(x)))]
        genre = self.city[person_index]
        
        kws = dict(mode='same', boundary=boundary)
        con = convolve2d(self.city == genre, self.kernel, **kws)
        neighs = convolve2d(self.city != 0, self.kernel, **kws)

        # new possible relocations are the ones that yet satisfy the "satisfaction" benchmark
        possible_relocations = (con / neighs >= 1 - self.T) & (self.city == 0)
        
        # list of index of suitable vacancies

        index_rc = [] 
        for x in np.arange(self.size):
            for y in np.arange(self.size):
                if possible_relocations[x, y]:
                    index_rc.append((x, y))

        if len(index_rc) == 0:
            # print("reached frozen configuration: no suitable relocation present")
            return None

        replace = index_rc[np.random.choice(np.arange(len(index_rc)))]
        self.city[person_index[0], person_index[1]] = 0
        self.city[replace] = genre

    def plot(self):
        """ This method plots the Schelling city  (in our case 2d system) visualizing it trought colors. """

        plt.style.use("ggplot")
        plt.figure(figsize=(12, 6))
        cmap = ListedColormap(['blue', 'white', 'red'])
        # blue is for -1, red for 1, ...
        plt.subplot(121)
        plt.axis('off')
        plt.title(f"size: {self.size} | t:{self.T} | empty:{self.p} ")
        # return plt.pcolor(self.city, cmap=cmap, edgecolors='w', linewidths=1)
        return plt.imshow(self.city, cmap=cmap, )
    
    

    """def sgregation_coefficient(self):
       # Returns the average segregation coefficient of the system, computed according to the (2)
       # formula in Gauvin et al paper
        
        s = 0
        coeff = 2 / ((self.size * self.size * (1 - self.p)) ** 2)
        for x in self.find_cluster_sizes():
            for n_c in x:
                s += n_c ** 2

        s = s * coeff
        return s"""

    def density_unwanted(self, boundary="fill", fillvalue = 0):
        """ This returns the two values of the densities of the unwanted locations for both
        types of agents. """

        kws = dict(mode='same', boundary=boundary, fillvalue = fillvalue) 
         
        # periodicity of boundaries is included or excluded thanks to our convolution
        neighs = convolve2d((self.city != 0).astype(int), self.kernel, **kws)
        con1 = convolve2d((self.city == 1).astype(int), self.kernel, **kws)
        # con_1 = convolve2d((self.city == -1).astype(int), self.kernel, **kws)

        densities = ((con1 / neighs < 1- self.T) & (self.city == 0))
        #print(densities)
        densities = densities.sum()
       # densities_1 = ((con_1 / neighs < 1-self.T) & (self.city == 0) ).sum()
        
        #densities of unwanted locations for agent +1 and agent -1
        return densities / self.vacancies

    def get_Ebc(self):

        if self.city.shape[0] < 3 or self.city.shape[1] < 3:
            raise ValueError("Input matrix is too small to compute periodic Moore neighbors sum")
            
        #extended_matrix = np.pad(self.city, ((1, 1), (1, 1)), mode='wrap') - if use this add mode = "valid"
        convolved_matrix = convolve2d(self.city, self.kernel, mode='same')
        result = np.multiply(self.city, convolved_matrix)
        return -result.sum()

    def get_Es(self):

        if self.city.shape[0] < 3 or self.city.shape[1] < 3:
            raise ValueError("Input matrix is too small to compute periodic Moore neighbors sum")
            
        #extended_matrix = np.pad(self.city, ((1, 1), (1, 1)), mode='wrap') - look comment get_ebc
        convolved_matrix = convolve2d(self.city, self.kernel, mode='same')
        convolved_matrix_squared = convolve2d(np.square(self.city), self.kernel, mode='same')
        result = -np.multiply(self.city, convolved_matrix) - self.T * (
            np.multiply(np.square(self.city), convolved_matrix_squared))
        return result.sum()

    def get_Es_normalized(self):
        "Return the normalized Energy, normalized by 4L2(1-p)"
        return self.get_Es() / (4 * self.size ** 2) * (1 - self.p)

    def get_Ebc_normalized(self):
        "Return the normalized Energy, normalized by 4L2(1-p)"
        return self.get_Ebc() / (4 * self.size ** 2) * (1 - self.p)

    def segregation_coefficient(self,matrix= True):
        """
        Find nearest neighbor clusters of +1 and -1 in a square matrix and compute
        the sum of the square of the number of components in each cluster.

        Args:
            matrix (numpy.ndarray): Square matrix representing the system, where
                                    1 represents cells occupied by group A, -1 represents
                                    cells occupied by group B, and 0 represents empty cells.

        Returns:
            int: Sum of the square of the number of components in each cluster.
            
        """
        matrix = self.city
        coeff  = 2/(self.size *self.size*(1-self.p))**2
        # Find connected components for group A (1)
        labeled_a, num_clusters_a = label(matrix == 1)

        # Find connected components for group B (-1)
        labeled_b, num_clusters_b = label(matrix == -1)

        # Initialize sum variable
        cluster_sum = 0

        # Calculate the sum of the square of the component counts for each cluster in group A
        for i in range(1, num_clusters_a + 1):
            cluster_indices = np.where(labeled_a == i)
            cluster_size = len(cluster_indices[0])
            cluster_sum += cluster_size ** 2

        # Calculate the sum of the square of the component counts for each cluster in group B
        for i in range(1, num_clusters_b + 1):
            cluster_indices = np.where(labeled_b == i)
            cluster_size = len(cluster_indices[0])
            cluster_sum += cluster_size ** 2

        return cluster_sum * coeff


    def real_space_renormalization(self):
        """ Method used to renormalize the model in presence of  high density of vacancies (p).
        Renormalization described in appendix A, paper Gauvin et al.
        """

        matrix = self.city
        updated_matrix = np.zeros((self.size, self.size))

        if self.size % 2 != 0:
            n = self.size - 1
        else:
            n = self.size

        for i in range(0, n, 2):
            for j in range(0, n, 2):
                neighbors = matrix[i:i + 2, j: j + 2]
                if (neighbors == 1).sum() > 2:
                    updated_matrix[i:i + 2, j: j + 2] = 1
                elif (neighbors == -1).sum() > 2:
                    updated_matrix[i:i + 2, j: j + 2] = -1
                elif (neighbors == 0).sum() > 2:
                    updated_matrix[i:i + 2, j: j + 2] = 0
                else:
                    updated_matrix[i:i + 2, j: j + 2] = matrix[i + 1][j + 1]
                # print(updated_matrix)

        self.city = updated_matrix
        
        
        return None 


    def renormalize(self):
        """
        Given a matrix of agents +1 , -1, and empty locations 0 renormalizes the matrix as follows:
        - if this site and its neighborhood comprise a majority of blue (resp. red) agents, the 2 × 2 square is replaced
        by a (single) blue (resp. red) agent;
        - if this site and its neighborhood consist of a majority of vacancies, the 2×2 square is replaced by a vacancy;
        - if there is no majority, the 4-site square is replaced by an agent of the same type as the bottom right site
        (or by a vacancy if that site is empty).
        """
        new_matrix = np.zeros_like(self.city)
        nrows = self.size
        ncols = self.size
        

        for i in range(nrows-1):
            for j in range(ncols-1):
                current = self.city[i:i+2, j:j+2]
                n_empty = np.sum(current == 0)
                n_blue = np.sum(current == 1)
                n_red = np.sum(current == -1)

                if n_empty >= 3:  # Majority of vacancies
                    new_matrix[i:i+2, j:j+2] = 0
                elif n_blue >= 3:  # Majority of blue
                    new_matrix[i:i+2, j:j+2] = 1
                elif n_red >= 3:  # Majority of red
                    new_matrix[i:i+2, j:j+2] = -1
                else:  # No majority
                    if current[1, 1] == 0:
                        new_matrix[i:i+2, j:j+2] = 0
                    elif current[1, 1] == 1:
                        new_matrix[i:i+2, j:j+2] = 1
                    else:
                        new_matrix[i:i+2, j:j+2] = -1
        
        
        self.city = new_matrix
        return None



