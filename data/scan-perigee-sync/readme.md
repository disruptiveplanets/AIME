## Aug 4

Created this scan today and staged it, but did not run it yet because the SC has too much going on. I'm going to wait until tomorrow, when gamma should be finished and much of vex and period will also be done.

## Aug 10

The code wasn't running fast enough so I killed it, tripled the CPU usage, and reran it. It should keep the chain and not kill anything, but the iteration numbering is wrong.

The result is that you should always use sample file number 1 (not zero) for analysis. This sped things up greatly.

## Aug 11

Pulled results so far (there are ten files which are still not completed) and pushed the directory to the SC. Then I ran the threshold extractors. I will need to rerun it for the remaining 10 later. I also made trim.py which sends all of the sample files to the right name.

These are the remaining 10.

peri-34.sh
peri-36.sh
peri-39.sh
peri-40.sh
peri-42.sh
peri-43.sh
peri-44.sh
peri-45.sh
peri-46.sh
peri-47.sh