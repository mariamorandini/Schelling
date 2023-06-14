import concurrent.futures

def run_dissatisfied_wrapper(model, i, j, k):
        model.run_dissatisfied()

def run_dissatisfied_multithreaded(T, P, N, C, Models):
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for i, t in enumerate(T):
            for j, p in enumerate(P):
                for n in N:
                    for k in C:
                        model = Models[i, j, k]
                        futures.append(executor.submit(run_dissatisfied_wrapper, model, i, j, k))

        # Wait for all threads to complete
        concurrent.futures.wait(futures)
