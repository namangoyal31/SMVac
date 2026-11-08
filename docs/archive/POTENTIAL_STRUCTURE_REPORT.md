# Potential Structure Report

## Parameters
- $M_h$ = 5 GeV
- $M_t$ = 105 GeV

## Analysis based strictly on $V(\phi)$ and $dV/d\phi$

### 1. Is $V(\phi)$ strictly monotonic?
**No**. The potential changes direction, as evidenced by $dV/d\phi$ changing sign.

### 2. Does $dV/d\phi$ change sign?
**Yes**. Sign changes observed at $\phi \approx$ ['2.54e+02', '1.08e+06'] GeV.

### 3. Does a barrier exist?
**Yes**. There is a region where $dV/d\phi$ changes from positive to negative, indicating a local maximum (barrier) before the potential turns downward.

# Potential Structure Report

## Parameters
- $M_h$ = 5 GeV
- $M_t$ = 105 GeV

## Analysis based strictly on $V(\phi)$ and $dV/d\phi$

### 1. Is $V(\phi)$ strictly monotonic?
**No**. The potential changes direction, as evidenced by $dV/d\phi$ changing sign.

### 2. Does $dV/d\phi$ change sign?
**Yes**. Sign changes observed at $\phi \approx$ ['2.54e+02', '1.08e+06'] GeV.

### 3. Does a barrier exist?
**Yes**. There is a region where $dV/d\phi$ changes from positive to negative, indicating a local maximum (barrier) before the potential turns downward.

### 4. Is there a second minimum?
**Yes**. There is a region where $dV/d\phi$ changes from negative to positive, indicating a local minimum at high field values.

## Plots
![Potential Plots](results/potential_structure_plot.png)

### 5. Stationary Points (Exact Analysis)
The exact stationary points were computed directly using root-finding algorithms (bisection) on $dV/d\phi = 0$. The potential incorporates the tree-level mass parameter $m^2 = \frac{1}{2} M_h^2$ to capture the low-field electroweak structure accurately.

The roots were verified independently by comparing the *numerical finite difference* of the potential with the *exact analytical derivative* (evaluating $\frac{d}{d\phi}\Delta\lambda$ inside the RGE using the $\beta$-functions chain-rule). Both methods agree within numerical precision.

* **Electroweak Minimum**: $\phi \approx 53.15$ GeV
  * $V(\phi_{EW}) \approx -6.48 \times 10^3 \text{ GeV}^4$
  * *(Note: At $M_h=5$ GeV, the extreme relative smallness of the tree-level quartic coupling compared to the quantum loop corrections shifts the EW minimum significantly away from $v \approx 246$ GeV).*
* **Barrier Maximum**: $\phi \approx 223.77$ GeV
  * $V(\phi_{bar}) \approx 1.96 \times 10^5 \text{ GeV}^4$
* **High-field (Second) Minimum**: $\phi \approx 1.08 \times 10^6$ GeV
  * $V(\phi_{hi}) \approx -5.93 \times 10^{19} \text{ GeV}^4$

#### Derived Quantities:
* **Barrier Height** ($\Delta V = V_{bar} - V_{EW}$): $2.02 \times 10^5 \text{ GeV}^4$
* **Second Minimum Depth** ($\Delta V = V_{hi} - V_{EW}$): $-5.93 \times 10^{19} \text{ GeV}^4$