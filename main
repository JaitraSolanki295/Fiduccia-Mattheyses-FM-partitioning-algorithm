import networkx as nx
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. INPUT DATA (Globals for easy function access)
# --------------------------------------------------
hyperedges = {'H1': {'a', 'b', 'c'}}
normal_edges = [('a', 'b'), ('c', 'd'), ('a', 'd'), ('a', 'e')]

# Initial partition for Pass 1
partA = {'a', 'b'}    # Blue partition
partB = {'c', 'd', 'e'}  # Green partition

# Custom positions
pos = {
    'a': (0, 0),
    'b': (0, 1),
    'c': (2, 1),
    'd': (2, 0),
    'e': (2, -1),
    'H1': (0.6, 0.5)
}
areas = {'a':2, 'b':4, 'c':1, 'd':4, 'e':5}

# --------------------------------------------------
# 2. HELPER FUNCTIONS
# --------------------------------------------------

def visualize_partition(hyperedges, normal_edges, partA, partB, pos=None, title="Partition Visualization"):
    """
    Visualize the given hypergraph + normal edges with partitions A and B.
    Partition A: Blue
    Partition B: Green
    Hyperedges: Orange
    """
    G = nx.Graph()

    for hname, nodes in hyperedges.items():
        G.add_node(hname, bipartite=1)
        for n in nodes:
            G.add_edge(hname, n)

    G.add_edges_from(normal_edges)

    if pos is None:
        pos = nx.spring_layout(G, seed=42)

    node_colors = []
    for n in G.nodes():
        if isinstance(n, str) and n.startswith('H'):
            node_colors.append('orange')
        elif n in partA:
            node_colors.append('#4A90E2')
        elif n in partB:
            node_colors.append('#50C878')
        else:
            node_colors.append('gray')

    plt.figure(figsize=(7, 5))
    nx.draw(G, pos, with_labels=True, node_color=node_colors,
            node_size=900, font_weight='bold', edge_color='gray')

    # if pos:
    #     plt.plot([1, 1], [-1.5, 1.5], 'k--', linewidth=2, label='Partition Boundary')

    plt.title(title, fontsize=13, fontweight='bold')
    plt.legend(loc='upper right')
    plt.axis('off')
    plt.show()

def calculate_fixed_balance_bounds(partA, partB, areas):
    """
    Calculates the FIXED balance criterion (BC) based on the initial partition.
    This BC will be applied to all subsequent passes as requested.
    """
    total_area = sum(areas.values())
    max_area = max(areas.values())
    area_A = sum(areas[n] for n in partA)
    area_B = sum(areas[n] for n in partB)

    # Calculate r based on the smaller area's ratio to total area
    smaller_area = min(area_A, area_B)
    r = smaller_area / total_area if total_area > 0 else 0.5

    # BC definition
    BC_low_raw = (r * total_area) - max_area
    BC_low = max(0, BC_low_raw)
    BC_high = (r * total_area) + max_area

    print("==================================================")
    print("FIXED BALANCE CRITERIA ESTABLISHED")
    print("==================================================")
    print(f"Total Area: {total_area}, Max Cell Area: {max_area}")
    print(f"Initial Area(A): {area_A}, Initial Area(B): {area_B}")
    print(f"Dynamic ratio r (FIXED for all passes): {r:.3f}")
    print(f"Balance Criterion Range (for Area A) (FIXED): [{BC_low:.2f}, {BC_high:.2f}]")
    print("==================================================\n")

    return BC_low, BC_high # Only return the fixed bounds

def compute_cell_gain(cell, hyperedges, normal_edges, partA, partB):
    """
    Compute gain Δg(c) = FS(c) - TE(c)
    """
    if cell in partA:
        own_part, other_part = partA, partB
    elif cell in partB:
        own_part, other_part = partB, partA
    else:
        return -float('inf'), 0, 0

    FS = TE = 0
    all_nets = [set(n) for n in hyperedges.values()] + [{u, v} for u, v in normal_edges]

    for net in all_nets:
        if cell not in net:
            continue
        others = net - {cell}
        if not others:
            continue

        # Case 1: all other cells in same partition -> uncut net -> contributes to TE
        if all(o in own_part for o in others):
            TE += 1
        # Case 2: all other cells in other partition -> cut net -> contributes to FS
        elif all(o in other_part for o in others):
            FS += 1

    gain = FS - TE
    return gain, FS, TE

def total_cutsize(nets, partA, partB):
    """
    Computes the total number of nets cut by the partition.
    """
    cut = 0
    for net in nets:
        inA = any(c in partA for c in net)
        inB = any(c in partB for c in net)
        if inA and inB:
            cut += 1
    return cut

# --------------------------------------------------
# 3. MAIN FM FUNCTION
# --------------------------------------------------

def fm_pass_detailed(hyperedges, normal_edges, partA, partB, areas, BC_low, BC_high, pos=None):
    """
    Performs a single, highly verbose pass of the FM algorithm,
    using the provided FIXED balance constraints.
    """
    print("STARTING FM PASS...")

    # --- 1. Initialization ---
    all_nets = [set(n) for n in hyperedges.values()] + [{u, v} for u, v in normal_edges]

    # Visualize the initial state of this specific pass
    visualize_partition(hyperedges, normal_edges, partA, partB, pos, title=f"Initial Partition (Start of FM Pass)")

    # Use the provided fixed BC for this pass
    print("==================================================")
    print("CURRENT PASS BALANCE CRITERIA (FIXED)")
    print(f"Balance Criterion Range (for Area A): [{BC_low:.2f}, {BC_high:.2f}]")
    print("==================================================\n")

    total_area = sum(areas.values())
    target_area = total_area / 2.0 # Used for tie-breaking

    current_partA = partA.copy()
    current_partB = partB.copy()
    all_cells = current_partA | current_partB
    locked_cells = set()
    free_cells = all_cells.copy()

    current_area_A = sum(areas[n] for n in current_partA)
    initial_cutsize = total_cutsize(all_nets, current_partA, current_partB)
    print(f"Initial Cutsize: {initial_cutsize}\n")

    move_history = []

    # Initialize cumulative gain tracking variables
    cumulative_gain = 0
    max_cumulative_gain = 0

    # History format: (Move Index, Cutsize, Cumulative Gain, Area A, Partition State)
    move_history.append((0, initial_cutsize, cumulative_gain, current_area_A, (partA.copy(), partB.copy())))

    # Best partition selection tracker
    best_partitions = [(0, current_area_A)] # Stores (Index, Area A) for iterations with Max Gain
    best_partition_index = 0

    # --- 2. Initial Gain Calculation (Sorted Output) ---
    print("Calculating initial gains for all cells...")
    gains = {}

    # Ensures output is in alphabetical order (a, b, c, d, e)
    sorted_free_cells = sorted(list(free_cells))
    for cell in sorted_free_cells:
        gain, fs, te = compute_cell_gain(cell, hyperedges, normal_edges, current_partA, current_partB)
        gains[cell] = gain
        print(f"Cell {cell}: FS = {fs}, TE = {te}, Gain = {gain:+.0f}")
    print("--------------------------------------------------")

    # --- 3. Iterative Moving ---
    num_moves_to_make = len(all_cells)
    for i in range(num_moves_to_make):
        print(f"\n----------------- MOVE {i+1} / {num_moves_to_make} -----------------")

        # --- a. Find best valid move ---
        # Sort free cells by gain (highest gain first)
        sorted_cells = sorted(
            [(c, g) for c, g in gains.items() if c in free_cells],
            key=lambda item: item[1],
            reverse=True
        )

        if not sorted_cells:
            print("No free cells left to move.")
            break

        cell_to_move = None
        move_gain = -float('inf')
        new_area_A = current_area_A # temp value

        for candidate_cell, candidate_gain in sorted_cells:
            print(f"Selected base cell -> {candidate_cell} (Gain = {candidate_gain:+.0f})")
            print(f"Checking BC before moving...")

            if candidate_cell in current_partA:
                new_area_A = current_area_A - areas[candidate_cell]
            else: # candidate_cell is in current_partB
                new_area_A = current_area_A + areas[candidate_cell]

            print(f"Area(A) = {new_area_A:.2f}")

            # FM only allows the move if it satisfies the FIXED BC
            if BC_low <= new_area_A <= BC_high:
                print(f"Balanced: True")
                print("--------------------------------------------------")
                print(f"BC Met -> Moving cell {candidate_cell}")
                cell_to_move = candidate_cell
                move_gain = candidate_gain
                break # Found the best valid move
            else:
                print(f"Balanced: False (Next cell will be checked)")
                # NOTE: If a cell violates BC, it is locked out of its partition for the pass.
                # Here, we only skip the move and proceed to the next highest gain cell.

        if cell_to_move is None:
            print(f"\nNo valid balanced move found. Stopping pass early.")
            break

        # --- b. Perform the move ---
        free_cells.remove(cell_to_move)
        locked_cells.add(cell_to_move)
        print(f"Locked cells: {sorted(list(locked_cells))}") # Sorted for consistent output

        if cell_to_move in current_partA:
            current_partA.remove(cell_to_move)
            current_partB.add(cell_to_move)
        else:
            current_partB.remove(cell_to_move)
            current_partA.add(cell_to_move)

        # Commit the new area
        current_area_A = new_area_A

        # Calculate CUMULATIVE GAIN
        cumulative_gain += move_gain

        # --- c. Record the new state and check for max gain ---
        current_cutsize = total_cutsize(all_nets, current_partA, current_partB)

        # Store comprehensive state in history
        move_history.append((i + 1, current_cutsize, cumulative_gain, current_area_A, (current_partA.copy(), current_partB.copy())))

        # Check if this move results in a NEW MAXIMUM CUMULATIVE GAIN or is a TIE
        if cumulative_gain > max_cumulative_gain:
            max_cumulative_gain = cumulative_gain
            best_partitions = [(i + 1, current_area_A)] # Reset list with new max
            print(f"New maximum cumulative gain found! (G_max = {max_cumulative_gain:+.0f})")
        elif cumulative_gain == max_cumulative_gain:
            best_partitions.append((i + 1, current_area_A)) # Add to tie list
            print(f"Tie in maximum cumulative gain found! (G_max = {max_cumulative_gain:+.0f})")

        # Update print statement with CUMULATIVE GAIN
        print(f"New State: Cutsize = {current_cutsize}, Cumulative Gain = {cumulative_gain:+.0f}, Area(A) = {current_area_A:.2f}")

        # --- d. Visualize ---
        visualize_partition(hyperedges, normal_edges, current_partA, current_partB, pos,
                            title=f"After Move {i+1}: '{cell_to_move}' (Cutsize: {current_cutsize}, Cum. Gain: {cumulative_gain:+.0f})")

        # --- e. Recompute gains for next iteration (Sorted Output) ---
        if free_cells:
            print(f"Recomputing gains for remaining unlocked cells...")
            gains.clear()

            # Ensures output is in alphabetical order (a, b, c, d, e)
            sorted_free_cells = sorted(list(free_cells))
            for cell in sorted_free_cells:
                # Gains must be computed based on the new partitions
                gain, fs, te = compute_cell_gain(cell, hyperedges, normal_edges, current_partA, current_partB)
                gains[cell] = gain
                print(f"Cell {cell}: FS = {fs}, TE = {te}, Gain = {gain:+.0f}")
            print("--------------------------------------------------")
        else:
            print("\nAll cells are locked.")

    # --- 4. Finalization: Apply Tie-Breaking Rule ---

    print(f"\n==================================================")
    print(f"PASS COMPLETE - FINAL SELECTION")
    print(f"==================================================")
    print(f"Maximum Cumulative Gain (G_max): {max_cumulative_gain:+.0f}")

    if max_cumulative_gain <= 0:
        # If no positive gain, return initial state (index 0)
        best_partition_index = 0
        print("No positive gain was achieved. Returning initial partition.")
    elif len(best_partitions) == 1:
        # No tie, simply use the index found
        best_partition_index = best_partitions[0][0]
        print(f"Unique Maximum Gain found at move {best_partition_index}.")
    else:
        # Tie-Breaking Logic: Select the partition closest to the ideal target area (best balance)
        print(f"Maximum Gain tie between moves: {[idx for idx, _ in best_partitions]}")

        min_balance_deviation = float('inf')
        final_selection_area = -1

        for idx, area_A in best_partitions:
            deviation = abs(area_A - target_area)
            print(f"Move {idx}: Area(A) = {area_A:.2f}, Deviation from Target ({target_area:.2f}) = {deviation:.2f}")

            if deviation < min_balance_deviation:
                min_balance_deviation = deviation
                best_partition_index = idx
                final_selection_area = area_A
            # Secondary tie-breaker: if deviations are equal, choose the earlier move
            elif deviation == min_balance_deviation and idx < best_partition_index:
                # This should not happen if the list is processed in order,
                # but included for robustness if needed.
                pass

        print(f"Tie broken: Selected move {best_partition_index} due to best balance (Area(A)={final_selection_area:.2f}).")


    # Retrieve the best partition state from history
    _ , final_cutsize, final_cumulative_gain, _, (final_partA, final_partB) = move_history[best_partition_index]

    print(f"\nResulting State after Pass Acceptance:")
    print(f"Initial cutsize was: {initial_cutsize}")
    print(f"Accepted Partition Index: {best_partition_index}")
    print(f"Final accepted cutsize: {final_cutsize}")

    # Return the best partition found, its cutsize, and the maximum gain achieved
    return final_partA, final_partB, final_cutsize, max_cumulative_gain


# --------------------------------------------------
# 4. EXECUTION WRAPPER
# --------------------------------------------------

def main():
    """Runs the full FM pass example using the predefined input data with iterative passes."""

    # --- CALCULATE FIXED BALANCE CRITERIA (Run once based on the initial partition) ---
    BC_low_fixed, BC_high_fixed = calculate_fixed_balance_bounds(partA, partB, areas)

    # --- ITERATIVE FM EXECUTION ---

    # Initial state variables
    current_A = partA.copy()
    current_B = partB.copy()
    pass_number = 1

    # Calculate initial cutsize
    all_nets = [set(n) for n in hyperedges.values()] + [{u, v} for u, v in normal_edges]
    current_cutsize = total_cutsize(all_nets, current_A, current_B)

    print("\n--- Initial Starting Partition ---")
    print(f"Partition A: {current_A}")
    print(f"Partition B: {current_B}")
    print(f"Initial Cutsize: {current_cutsize}\n")

    while True:
        print("\n\n##################################################")
        print(f"                 STARTING PASS {pass_number}")
        print("##################################################")

        # Run the detailed FM pass. The starting partition is the result of the last accepted pass.
        # Note: BC_low_fixed and BC_high_fixed are used in every pass.
        next_A, next_B, next_cutsize, max_gain = fm_pass_detailed(
            hyperedges, normal_edges, current_A, current_B, areas, BC_low_fixed, BC_high_fixed, pos
        )

        print("\n--- Pass Results Summary ---")
        print(f"Maximum Cumulative Gain in Pass {pass_number}: {max_gain:+.0f}")
        print(f"Accepted Partition A: {next_A}")
        print(f"Accepted Partition B: {next_B}")
        print(f"Accepted Cutsize: {next_cutsize}")

        # The stopping condition: If the maximum cumulative gain is 0 or negative,
        # it means no move sequence resulted in an improvement.
        if max_gain <= 0:
            print("\n==================================================")
            print(f"STOPPING CRITERIA MET: Maximum gain in Pass {pass_number} was {max_gain:+.0f}.")
            print("No further improvement is possible. Final partition determined.")
            print("==================================================")
            break

        # If max_gain > 0, the partition was improved. Update state for the next pass.
        print(f"\n--- Improvement found. Proceeding to Pass {pass_number + 1} ---")
        current_A = next_A
        current_B = next_B
        current_cutsize = next_cutsize
        pass_number += 1

    # Final Output after the loop breaks
    print("\n\n##################################################")
    print("           FINAL STABLE PARTITION")
    print("##################################################")
    print(f"Total FM Passes Executed: {pass_number}")
    print(f"Final Partition A: {current_A}")
    print(f"Final Partition B: {current_B}")
    print(f"Final Cutsize: {current_cutsize}")

    # Visualize the FINAL stable partition one last time
    visualize_partition(hyperedges, normal_edges, current_A, current_B, pos,
                        title=f"FINAL STABLE PARTITION (Cutsize: {current_cutsize})")


if __name__ == "__main__":
    main()      
