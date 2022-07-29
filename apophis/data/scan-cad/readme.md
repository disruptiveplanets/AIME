## April 26

Started with the tumbling asteroid with a=263 (sqrt(3/5) times [340](https://echo.jpl.nasa.gov/asteroids/brozovic.etal.apophis.2018.pdf)) and the log normal rho uncertainty. I've used the DRC as e-3 on both sides for testing's sake, but plan to extend it.

## April 28

The code was taking longer than I expected, so I copied one and named it cad-99 to pull it. It was highly degenerate, so I went back into the code to figure out why. I found a bug having to do with the initial orientation, so I reran cad-99 which is a duplicate of cad-00, but with the correct orientation and with the parameters switched from (gamma, phi) to (gamma+phi, gamma-phi).

## April 29

Had to kill the fits to run other things, but it's clear that degeneracy is a problem. Removed the + and - from main.py to avoid mistakes.

## May 6

Pushed today

## May 9

Pulled today. For some reason, lots of them failed. So I'm re-pushing those. Of those that succeeded, there are still some issues. So I'm working on that.