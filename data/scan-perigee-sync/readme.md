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

Then I realized that some had errors:

peri-10.sh
peri-11.sh
per-12.sh
peri-21.sh
peri-23.sh
peri-29.sh
peri-31.sh
peri-37.sh
peri-41.sh

and I reran those too. However, they didn't look too bad. The only runs I was sus of were

11 (weak data)
20 (weak data)
22 (no data)
25 (bad burnin)
27 (bad burnin)
28 (bad burnin)
29 (bad burnin)
30 (bad burnin)
31 (weak data)
32 (weak data)
35 (weak data)
37 (no data)
38 (weak data)
40 (no data)


given the PNGs. So I also reran

20
22
25
27
28
30
32
35
38

All of them were on reload, except 37, which had no data. So total running are 

peri-10.sh
peri-11.sh
peri-12.sh
peri-20.sh
peri-21.sh
peri-22.sh
peri-23.sh
peri-25.sh
peri-27.sh
peri-28.sh
peri-29.sh
peri-30.sh
peri-31.sh
peri-32.sh
peri-34.sh
peri-35.sh
peri-36.sh
peri-37.sh
peri-38.sh
peri-39.sh
peri-40.sh
peri-41.sh
peri-42.sh
peri-43.sh
peri-44.sh
peri-45.sh
peri-46.sh
peri-47.sh

Tomorrow I will 

* Check for errors in the log files
* Ask how many have converged and how quickly
* Redo thresholds for them if any have passed
* Edit paper

## Aug 12

No errors. Thresholds are not yet finished. Here's what's still running:

peri-23.sh
peri-29.sh                  
peri-31.sh                  
peri-34.sh                  
peri-36.sh                  
peri-37.sh                  
peri-40.sh                  
peri-41.sh                  
peri-42.sh                  
peri-43.sh                  
peri-44.sh                  
peri-46.sh                  
peri-47.sh

which means these finished:

peri-10.sh (good)
peri-11.sh (good)
peri-12.sh (good)
peri-20.sh (fine)
peri-21.sh (good)
peri-22.sh (good)
peri-25.sh (good)
peri-27.sh (good)
peri-28.sh (good)
peri-30.sh (good)
peri-32.sh (good)
peri-35.sh (fine)
peri-38.sh (weird degeneracy)
peri-39.sh (fine)
peri-45.sh (fine)

I re-terminated them and pulled their data to confirm that they were acceptable. But in so doing, I accidentally killed all the other processes and had to reload them. However, this did not induce any errors. The number of iterations already completed for all processes were between 7,000 and 30,000. I will manually kill them at 100,000.

To check the completed files, I deleted the directories which are still running and extracted plots. They look good, but I had to change 33 to adjust the burnin manually and saved it as peri-33-1.