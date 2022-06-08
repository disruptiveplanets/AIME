## Jun 2

I started extraction of uncertainties as a function of observational precision. Division of 49 m.

## June 3

I killed the process and switched to 99 m. It was taking too long. (2 hours per task). I'm running all scans now.

## June 6

I edited the code so that it would not do reruns of data already generated, and would not kill itself when it hit a failed run. Pulled some of the scans.

I also realized that the output is pure uncertainty, not uncertainty ratio. So I changed it to uncertainty ratio.

After my meeting with Julien, I changed the prior to 0.25-3 and adjusted the number of finite elements used. Then I ran all of these on the SC.

## June 8

I finished the DOF scan and found that fewer DOFs leads to less uncertainty, not more. Now I'm rerunning with only a few DOFs, averaging over many runs, and booting all those that fail. This is being parallelized.