* **Run 0:**
```
CADENCE = 3600
spin = [0.00012, 0.00012, 0.00012]
jlms = [5.972e24, 5.972e22, -5.972e22, 4.972e22]
theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)
theta_start = (
    0.1, 1.0e6, 1.0e5, -5.0e5,
)
theta_range = (
    (0, 2 * np.pi), (0.5e6, 2e6), (0.5e5, 2.0e5), (-2.5e5, -1.0e6),
)
```
* **Run 1:**
```
CADENCE = 3600
spin = [0.00012, 0.00022, 0.00032]
jlms = [5.972e24, 5.972e22, -5.972e22, 4.972e22]
theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)
theta_start = (
    0.1, 1.0e6, 1.0e5, -5.0e5,
)
theta_range = (
    (-np.pi, -np.pi), (0.5e6, 2e6), (0.5e5, 2.0e5), (-2.5e5, -1.0e6),
)
```
Also randomized data, `REL_SIGMA=0.1`
* **Run 2:**
```
CADENCE = 120
spin = [0.00012, 0.00022, 0.00032]
jlms = [5.972e24, 0, 0, 4.972e22]
theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)
theta_start = (
    0.1, 1.0e6, 1.0e5, -5.0e5,
)
theta_range = (
    (-np.pi, np.pi), (0.5e6, 2e6), (0.5e5, 2.0e5), (-2.5e5, -1.0e6),
)
```
Also randomized data, `REL_SIGMA=0.1`
