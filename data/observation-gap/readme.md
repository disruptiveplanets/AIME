## Jan 27

Created today. I pushed things to the supercomputer as well.

## Jan 29

Completed today. I pulled things from the supercomputer after making a slight adjustment to display to make sure it could handle gaps. Everything worked fine, but I got some slightly unexpected results. So I changed the gap range from (0, 1) hr to (0, 2) hr and I upped the amount of data that was accepted (before, it was too low). Now it should yield some more interesting data.

## Jan 31

Completed today. Pulled everything but gap 14. Things still don't look quite right, so I extended to 3 hours and tried again.

## Feb 1

Reran the program because of a mistake I made (forgot to change the mask size)

## Feb 2

I pulled the results today, but the data still doesn't look right. It has these weird tiered increase structures, that sometimes line up with integer hours. So I ran the benchmark test with the same system to see if it shows the same features. It does not. I re-pulled the true test and got the same thing. So I think it's real. I also fixed the summarizing figure, with the bands indicating divergence in the spin vectors. I wrote up the observation gap section as well.