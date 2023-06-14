import numpy as np

def optimize_np(T, P, N, C, Models):
    """ in this version of the optimized method we store also all the intermediate stages 
    of evolution of the class instances """

    T_size = len(T)
    P_size = len(P)
    C_size = len(C)
    N_size = len(N)

    # Create an empty array to store the city values
    Models_storage = np.empty((T_size, P_size, C_size, N_size), dtype=object)

    for i in range(T_size):
        for j in range(P_size):
            for n in range(N_size):
                for k in range(C_size):
                    model = Models[i, j, k]
                    model.run_dissatisfied()
                    Models_storage[i, j, k, n] = model.city

    return Models_storage