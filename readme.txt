DN52esmc:subread ziliangq$ /usr/bin/time -lp ./a.out ../nm_publication/chr1.fa ../nm_publication/SRR1197524.1M.fastq
current time in ms  1429298269
TOTAL:165462759 MULTI:83787846
current time in ms  1429298440
current time in ms  1429298536 (1m36s! for a 1Mx200bp reads against whole hg, 200M/100s=2MBps, single core, 3 times compared to original subread)
done with hits:0
real       267.45
user       255.44
sys          8.75
7138213888  maximum resident set size
         0  average shared memory size
         0  average unshared data size
         0  average unshared stack size
   3063600  page reclaims
         0  page faults
         0  swaps
         1  block input operations
         1  block output operations
         0  messages sent
         0  messages received
         0  signals received
         8  voluntary context switches
    160332  involuntary context switches

Usaully, we have 20M reads x 200bp, it would take 40min/sample/thread. We need to make the program 40 times faster! 
Strategy, 
1, remove highly genomic repeats during mapping
2, we usually have 2x2 CPU cores
it is still possible to make 1min computing




