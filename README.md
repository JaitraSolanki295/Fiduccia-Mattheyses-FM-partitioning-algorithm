# FM algorithm
## Gain Function

In the FM algorithm, the gain of each movable cell is calculated to determine whether moving the cell to the opposite partition will reduce the cutsize.

The gain of a cell **c** is defined as

\[
Δg(c)=FS(c)-TE(c)
\]

where

- **FS(c)** = Number of cut nets connected to **c** that become uncut after moving the cell.
- **TE(c)** = Number of uncut nets connected to **c** that become cut after moving the cell.

A positive gain indicates that moving the cell will reduce the cutsize and is therefore preferred.


The cell having the **maximum gain** is selected for movement, provided it satisfies the balance criterion.

---

## Ratio Factor

The **Ratio Factor (r)** represents the relative balance between the two partitions based on their total cell area.

\[
   r= (area(A))/(area(A)+area(B))
\]

where

- **Area(A)** = Total area of cells in Partition A
- **Area(B)** = Total area of cells in Partition B

The ratio factor is used to prevent all cells from clustering into a single partition and helps maintain balanced partition sizes.

---

### Balance Criterion (BC)

The FM algorithm maintains balanced partitions by ensuring that every cell movement satisfies the balance criterion.

The allowable range for the area of Partition A is

\[
**BC** =>   [r * T - Amax] < **area [A]** <   [r * T + Amax]
\]

where,

- **T** = Total area of all cells
- **Amax** = Maximum cell area

Only the cell movements satisfying this condition are accepted.

---

# FM Algorithm Steps

### Step 1 : Compute the Balance Criterion (BC)

Compute the balance criterion using the ratio factor and the area of all cells.

---

### Step 2 : Compute Cell Gain

Calculate the gain of every movable cell using

\[
Δg=FS(c)-TE(c)
\]

---

### Step 3 : Select the Base Cell

Choose the unlocked cell having the **maximum gain** that also satisfies the **Balance Criterion (BC)**.

Move the selected cell to the opposite partition.

---

### Step 4 : Lock the Base Cell

Lock (fix) the selected cell so that it cannot be moved again during the current pass.

Recompute the gains of all remaining unlocked cells.

---

### Step 5 : Continue Until All Cells are Locked

Repeat the gain computation and cell movement until every cell has been locked.

---

### Step 6 : Determine the Best Partition

After all cells are locked,

- Calculate the cumulative gain after every move.
- Identify the partition corresponding to the **maximum cumulative gain**.
- Select this partition as the final partition for the current pass.

---

### Step 7 : Execute the Next Pass

If

\[
G_{max}>0
\]

then

- Unlock all fixed cells.
- Use the best partition obtained in the previous pass as the initial partition for the next pass.
- Repeat the complete FM algorithm.

Otherwise,

- Stop the algorithm since no further improvement is possible.

---




