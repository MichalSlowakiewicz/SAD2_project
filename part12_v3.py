# This code covers Task 1 and Task 2 from Part 1

# -*- coding: utf-8 -*-
from __future__ import print_function  # import for using standard 'print()' from Python 3 in Python 2
import networkx as nx
import boolean
import random
import itertools
import pandas as pd
import os

# Note: class BN bases on the one introduced during lab classes in 'bn_exercise.py' code
class BN(object):
    __bool_algebra = boolean.BooleanAlgebra()

    def __init__(self, list_of_nodes, list_of_functions):

        self.num_nodes = len(list_of_nodes)
        self.node_names = list_of_nodes

        self.list_of_nodes = []
        for node_name in list_of_nodes:
            node = self.__bool_algebra.Symbol(node_name)
            self.list_of_nodes.append(node)

        self.functions = []
        for fun in list_of_functions:
            self.functions.append(self.__bool_algebra.parse(fun).simplify())

    def __int_to_state(self, x):
        """
        Helper method for converting a non-negative integer into a state in the form of a tuple of 0s and 1s.
        Args:
            x (int): A state number
        Returns:
            tuple[int, ...]: A tuple of 0s and 1s representing the Boolean network state.
        """
        binary_str = format(x, '0' + str(self.num_nodes) + 'b')
        state = [int(char) for char in binary_str]
        return tuple(state)

    def evaluate_node(self, node_idx, state_tuple):
        """
        Calculates the boolean value of a certain node given the network.
        """
        substitution = {}
        for i in range(self.num_nodes):
            if state_tuple[i] == 1:
                substitution[self.list_of_nodes[i]] = self.__bool_algebra.TRUE
            else:
                substitution[self.list_of_nodes[i]] = self.__bool_algebra.FALSE

        result = self.functions[node_idx].subs(substitution).simplify()

        return 1 if result == self.__bool_algebra.TRUE else 0

    def get_neighbor_states(self, state):
        """
        Computes the states reachable from the given state in one step of asynchronous update.
        Args:
            state (tuple[int, ...]): A tuple of 0s and 1s representing the Boolean network state.
        Returns:
            set[tuple[int, ...]]: A set of tuples of 0s and 1s representing the Boolean network states reachable
                in one step from the given state.
        """
        neighbors = set()

        for i in range(self.num_nodes):
            new_val = self.evaluate_node(i, state)

            # creating new state
            next_state_list = list(state)

            if next_state_list[i] != new_val:
                next_state_list[i] = new_val
                neighbors.add(tuple(next_state_list))
            else:
                neighbors.add(state)

        return neighbors

    def generate_state_transition_system(self):
        """
        Generates the asynchronous state transition system of the Boolean network.
        Returns:
            nx.DiGraph: NetworkX DiGraph object representing the asynchronous state transition system.
        """
        G = nx.DiGraph()

        for i in range(2 ** self.num_nodes):
            current_state = self.__int_to_state(i)
            G.add_node(current_state)

            neighbors = self.get_neighbor_states(current_state)
            for neighbor in neighbors:
                G.add_edge(current_state, neighbor)
        return G

    def get_next_state_synchronous(self, state):
        """
        Computes the next state in synchronous mode.
        """
        next_state_list = []
        for i in range(self.num_nodes):
            next_state_list.append(self.evaluate_node(i, state))
        return tuple(next_state_list)

    def simulate_trajectory(self, start_state, steps, mode='sync', sampling_freq=1):
        """
        Simulates single trajectory for dataset.
        """
        trajectory = []
        current_state = start_state

        trajectory.append(current_state)


        for t in range(1, steps * sampling_freq + 1):
            if mode == 'sync':
                current_state = self.get_next_state_synchronous(current_state)
            else:  # asynchronous mode
                # drawing the index of the node
                node_to_update = random.randint(0, self.num_nodes - 1)
                new_val = self.evaluate_node(node_to_update, current_state)
                next_s = list(current_state)
                next_s[node_to_update] = new_val
                current_state = tuple(next_s)

            if t % sampling_freq == 0:
                trajectory.append(current_state)

        return trajectory

    def analyze_attractor_structure(self, start_state, max_steps=10000):
        """
        Simulates the network and returns two objects: (transient, attractor).
        'transient' represent transient states and 'attractor' represent states belonging to attractor.
        Simulation is conducted by making 10k steps in the network, which is sufficient taking into account that 16 is the maximal number of nodes.
        """
        history = []
        seen = {}
        curr = start_state


        for i in range(max_steps):
            if curr in seen:
                # found cycle
                start_index = seen[curr]
                transient_part = history[:start_index]
                attractor_part = history[start_index:]
                return transient_part, attractor_part

            seen[curr] = i
            history.append(curr)
            curr = self.get_next_state_synchronous(curr)

        # if no cycles were found
        return history, []

    def simulate_trajectory_async(self, start_state, steps, sampling_freq=1):
        """
        Generates dataset for Task 2 (asynchronous).
        """
        trajectory = [start_state]
        current_state = start_state

        # simulation loop
        for t in range(1, steps * sampling_freq + 1):
            # drawing node for update
            node_to_update = random.randint(0, self.num_nodes - 1)
            new_val = self.evaluate_node(node_to_update, current_state)

            # creating new state
            next_s = list(current_state)
            next_s[node_to_update] = new_val
            current_state = tuple(next_s)

            # saving results once every 'sampling_freq' step
            if t % sampling_freq == 0:
                trajectory.append(current_state)

        return trajectory

def generate_random_bn_strings(num_nodes, max_parents=3):
    """
    Generates nodes' names and random functions in text format.
    """

    node_names = ['x{}'.format(i + 1) for i in range(num_nodes)]
    functions = []

    for i in range(num_nodes):
        k = random.randint(1, min(max_parents, num_nodes))
        parents = random.sample(node_names, k)

        parent_combinations = list(itertools.product([0, 1], repeat=k))
        true_combinations = [c for c in parent_combinations if random.random() > 0.5]

        if not true_combinations:
            functions.append("0")
        elif len(true_combinations) == len(parent_combinations):
            functions.append("1")
        else:
            clauses = []
            for combo in true_combinations:
                literals = []
                for parent_idx, val in enumerate(combo):
                    p_name = parents[parent_idx]

                    literals.append(p_name if val == 1 else "~{}".format(p_name))
                clauses.append("(" + " & ".join(literals) + ")")

            full_func = " | ".join(clauses)
            functions.append(full_func)

    return node_names, functions

def save_to_file(data, colnames, filename):
    """
    Saves dataset in format compliant with BNFinder2.
    """
    if not data: return

    # creating dataframe
    df = pd.DataFrame(data, columns=colnames)

    df_transposed = df.T

    # name for column with variable names
    df_transposed.index.name = "GeneName"

    # saving dataframe as csv/txt
    df_transposed.to_csv(filename, sep='\t', header=True, index=True)

# main pipeline
if __name__ == "__main__":

    # output folder for datasets
    output_dir = "Project_Data_Final_5"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # chosen sizes of boolean networks
    sizes_to_test = [5, 7, 10, 16]
    # numbers of starts for each network
    NUM_STARTS = 300

    # loop resposnsible for generating datasets
    for N in sizes_to_test:

        # number of steps for trajectories with asynchronous updates
        ASYNC_STEPS = N * 5

        # generating network
        names, funcs = generate_random_bn_strings(num_nodes=N, max_parents=3)
        bn = BN(names, funcs)

        # saving ground truth
        with open("{}/network_def_{}.txt".format(output_dir, N), "w") as f:
            for n_name, f_func in zip(names, funcs):
                f.write("{} = {}\n".format(n_name, f_func))

        # ---------------------------------
        # --- SINGLE TRAJECTORY SECTION ---
        # ---------------------------------

        # in this section datasets containing only one trajectory are generated
        # for both synchronous and asynchronous mode we generate 4 datasets:
        # entire trajectory, only the beginning (transient + small part of attractor),
        # only the end (small part of transient + full attractor), entire trajectory sampled every 3rd step
        start_node_single = tuple([random.randint(0, 1) for _ in range(N)])

        # SYNCHRONOUS DATASETS
        # We split the simulation history into 'transient' and 'attractor' part
        # to create specialized datasets (with different proportion of transient and attractor states).
        s_trans, s_attr = bn.analyze_attractor_structure(start_node_single)
        if s_attr:
            s_full = s_trans + (s_attr * 5)
            save_to_file(s_full, names, "{}/data_size{}_sync_SINGLE_full.txt".format(output_dir, N))
            save_to_file(s_trans + s_attr[:3], names, "{}/data_size{}_sync_SINGLE_trans_heavy.txt".format(output_dir, N))
            save_to_file(s_trans[-3:] + s_attr, names, "{}/data_size{}_sync_SINGLE_attr_heavy.txt".format(output_dir, N))
            save_to_file(s_full[::3], names, "{}/data_size{}_sync_SINGLE_sampled_f3.txt".format(output_dir, N))

        # ASYNCHRONOUS DATASETS
        # Since the exact attractor detection in asynchronous setting is hard, we assume that
        # early steps are transient and late steps belong to attractor.
        s_raw_async = bn.simulate_trajectory_async(start_node_single, steps=ASYNC_STEPS)
        save_to_file(s_raw_async, names, "{}/data_size{}_async_SINGLE_full.txt".format(output_dir, N))
        save_to_file(s_raw_async[:int(ASYNC_STEPS*0.3)], names, "{}/data_size{}_async_SINGLE_trans_heavy.txt".format(output_dir, N))
        save_to_file(s_raw_async[-int(ASYNC_STEPS*0.3):], names, "{}/data_size{}_async_SINGLE_attr_heavy.txt".format(output_dir, N))
        save_to_file(s_raw_async[::3], names, "{}/data_size{}_async_SINGLE_sampled_f3.txt".format(output_dir, N))

        # -----------------------------------
        # --- MULTIPLE TRAJECTORY SECTION ---
        # -----------------------------------

        # in this section datasets containing only multiple (NUM_STARTS) trajectories are generated
        # for both synchronous and asynchronous mode we generate 4 datasets:
        # entire trajectory, only the beginning (transient + small part of attractor),
        # only the end (small part of transient + full attractor), entire trajectory sampled every 3rd step

        # SYNCHRONOUS DATASETS
        sync_full = []  # entire trajectory (transient + attractor)
        sync_trans_heavy = []  # only the beginning (transient + small part of attractor)
        sync_attr_heavy = []  # only the end (small part of transient + full attractor)
        sync_sampled_f3 = []  # entire trajectory sampled every 3rd step

        # drawing initial state
        current_start = tuple([random.randint(0, 1) for _ in range(N)])

        # 'NUM_STARTS' independent runs
        for _ in range(NUM_STARTS):
            trans, attr = bn.analyze_attractor_structure(current_start)

            # skipping if no cycle was found
            if not attr:
                current_start = tuple([random.randint(0, 1) for _ in range(N)])
                continue

            full_traj = trans + (attr * 5)
            sync_full.extend(full_traj)
            sync_trans_heavy.extend(trans + attr[:3])
            sync_attr_heavy.extend(trans[-3:] + attr)
            sync_sampled_f3.extend(full_traj[::3])

            # perturbation at the end of the trajectory:
            # we take the last state of the trajectory and change one gene before simulating the next trajectory
            last_state = list(full_traj[-1])
            idx_to_flip = random.randint(0, N - 1)
            last_state[idx_to_flip] = 1 - last_state[idx_to_flip]
            current_start = tuple(last_state)

        save_to_file(sync_full, names, "{}/data_size{}_sync_full.txt".format(output_dir, N))
        save_to_file(sync_trans_heavy, names, "{}/data_size{}_sync_transient_heavy.txt".format(output_dir, N))
        save_to_file(sync_attr_heavy, names, "{}/data_size{}_sync_attractor_heavy.txt".format(output_dir, N))
        save_to_file(sync_sampled_f3, names, "{}/data_size{}_sync_sampled_f3.txt".format(output_dir, N))

        # ASYNCHRONOUS DATASETS
        async_full = []  # entire trajectory
        async_trans_heavy = []  # first 50 steps
        async_attr_heavy = []  # last 50 steps
        async_sampled_f3 = []  # entire trajectory sampled every 3rd step

        # 'NUM_STARTS' independent runs
        for _ in range(NUM_STARTS):
            start_node = tuple([random.randint(0, 1) for _ in range(N)])

            # simulating asynchronous steps
            raw_traj = bn.simulate_trajectory_async(start_node, steps=ASYNC_STEPS, sampling_freq=1)

            # Dataset 1: taking full trajectory
            async_full.extend(raw_traj)

            # Dataset 2: sampling every third observation from trajectory
            sampled_part = raw_traj[::3]
            async_sampled_f3.extend(sampled_part)

            # Dataset 3: taking 50 first observations. It represents the immediate response to initial conditions - should contain a lot of transient.
            async_trans_heavy.extend(raw_traj[:int(ASYNC_STEPS*0.3)])

            # Dataset 4: taking 50 last observations. It represents settled system - it should contain settle in an attractor.
            async_attr_heavy.extend(raw_traj[-int(ASYNC_STEPS*0.3):])

        save_to_file(async_full, names, "{}/data_size{}_async_full.txt".format(output_dir, N))
        save_to_file(async_trans_heavy, names, "{}/data_size{}_async_trans_heavy.txt".format(output_dir, N))
        save_to_file(async_attr_heavy, names, "{}/data_size{}_async_attr_heavy.txt".format(output_dir, N))
        save_to_file(async_sampled_f3, names, "{}/data_size{}_async_sampled_f3.txt".format(output_dir, N))
