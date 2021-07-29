
### July 27, 2021
<table><tr>
<th width=800px>
Run 0:

```
CADENCE = 3600
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
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
</th><th>
<img src="run-0-20000.png" alt="corner plot" width="400"/>
</th></tr></table>

The posterior for psi is too clustered up near zero because the true value is at the edge of the prior. Need to change the prior. Also, I don't like my error analysis. I'm going to change it to be proportional to the value itself.



---
### July 28, 2021

I'm running a run with the new prior and a new spin to test the new uncertainty (1) and one with more cadence to test that effect, with only J10 to see if it makes a difference (2).

<table><tr>
<th width=800px>
Run 1:

```
CADENCE = 3600
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
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
REL_SIGMA=0.1
```
</th><th>
<img src="run-1-5000.png" alt="corner plot" width="400"/>
</th></tr><tr>
<th width=800px>
Run 2:

```
CADENCE = 120
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
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
REL_SIGMA=0.1
```
</th><th>
<img src="run-2-5000.png" alt="corner plot" width="400"/>
</th></tr></table>

The degeneracy remains in run 1 but is broken in run 2. Probably because the `jlm[1]` and `jlm[2]` have equal magnitude. I'll test that in the future. The other option is that the degeneracy is only broken if `jlm[3]` is big enough, or perhaps more cadence is needed to break degeneracy.

Neither corner plot is very good, probably due to lack of data. I could always run the fit for longer. In the future, I'll run for 7500 iterations, not 5000.

The fits for run 2 are nowhere near the correct value, so I'm going to change the error analysis again. This time it simulates the spin vector having an incorrect angular component, but a correct radial component. The error bars on large values are therefore smaller now. Hopefully that produces better fit behavior.




---
### July 29, 2021

I'm running the same as runs 1 and 2 but with the new error to detect improvement (3, 4). Then I'm running the same as run 0 for the same reason (5). To see if the degeneracy is broken by having `|jlm[1]| != |jlm[2]|`, I'm running run 0 with a new `jlm[2]` (6). Then the same as run 6 but with low cadence to test that effect in isolation.

<table><tr>
<th width=800px>
Run 3:

```
CADENCE = 3600
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
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
SIGMA=0.2
```
</th><th>
<img src="run-3-5000.png" alt="corner plot" width="400"/>
</th></tr><tr>
<th width=800px>
Run 4:

```
CADENCE = 120
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
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
SIGMA=0.2
```
</th><th>
<img src="run-4-5000.png" alt="corner plot" width="400"/>
</th></tr><tr>
<th width=800px>
Run 5:

```
CADENCE = 3600
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
spin = [0.00012, 0.00012, 0.00012]
jlms = [5.972e24, 5.972e22, -5.972e22, 4.972e22]
theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)
theta_start = (
    0.1, 1.0e6, 1.0e5, -5.0e5,
)
theta_range = (
    (0, 2 * np.pi, (0.5e6, 2e6), (0.5e5, 2.0e5), (-2.5e5, -1.0e6),
)
SIGMA=0.2
```
</th><th>
<img src="run-5-5000.png" alt="corner plot" width="400"/>
</th></tr><tr>
<th width=800px>
Run 6:

```
CADENCE = 3600
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
spin = [0.00012, 0.00012, 0.00012]
jlms = [5.972e24, 5.16143e22, -7.421e22, 4.972e22]
theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)
theta_start = (
    0.1, 1.0e6, 1.0e5, -5.0e5,
)
theta_range = (
    (0, 2 * np.pi, (0.5e6, 2e6), (0.5e5, 2.0e5), (-2.5e5, -1.0e6),
)
SIGMA=0.2
```
</th><th>
<img src="run-6-5000.png" alt="corner plot" width="400"/>
</th></tr><tr>
<th width=800px>
Run 7:

```
CADENCE = 120
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
spin = [0.00012, 0.00012, 0.00012]
jlms = [5.972e24, 5.16143e22, -7.421e22, 4.972e22]
theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)
theta_start = (
    0.1, 1.0e6, 1.0e5, -5.0e5,
)
theta_range = (
    (0, 2 * np.pi, (0.5e6, 2e6), (0.5e5, 2.0e5), (-2.5e5, -1.0e6),
)
SIGMA=0.2
```
</th><th>
<img src="run-6-5000.png" alt="corner plot" width="400"/>
</th></tr></table>
