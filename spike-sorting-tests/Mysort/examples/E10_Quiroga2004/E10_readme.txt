In his 2004 publication "Unsupervised spike detection and sorting with 
wavelets and superparamagnetic clustering" Quiroga et al published a 
benchmark data set for spike sorting. The data set consists of 4 sub data
sets (Easy 1 and 2 and Difficult 1 and 2). In this example we show, how
this dataset can be analysed and the available ground truth information
(it is simulated data, so we know where the spikes are) used to compute
the performance of the sorting.

You have to download the benchmark data that comes with the spike sorter 
"wave_clus". Then, set the path to the folder containing the files in 
the "E10_quirogaDataPath.m" file accordingly.

Run E10_main.m and you should see the following output:

    Evaluation reported in 2004 publication Quiroga et al:
    Easy 1    |#Spks|#Ovps|#FN|#FNOv|#FP|clas.Spks|clas.Err |
    ----------|---------------------------------------------|
    Noise 0.05| 3514   785  17   193 711      2729        1 |
    Noise 0.10| 3522   769   2   177  57      2753        5 |
    Noise 0.15| 3477   784 145   215  14      2693        5 |
    Noise 0.20| 3474   796 714   275  10      2678       12 |
    Noise 0.25|   -1    -1  -1    -1  -1      2586       64 |
    Noise 0.30|   -1    -1  -1    -1  -1      2629      276 |
    Noise 0.35|   -1    -1  -1    -1  -1      2702      483 |
    Noise 0.40|   -1    -1  -1    -1  -1      2645      741 |
    ----------|---------------------------------------------|
    Total     |13983  3130 874   856 788     21415     1587 |
    Preprocessing ... C_Easy1_noise02.mat
    ...preprocessing done.
    Warning, noise covariance matrix ill conditioned. Applying subspace loading. (cond = 15749.0898)
    Processing Chunk  1/15 (  7%)
    Processing Chunk  2/15 ( 13%, est:    0 min:04 sek):
    Processing Chunk  3/15 ( 20%, est:    0 min:04 sek):
    Processing Chunk  4/15 ( 27%, est:    0 min:03 sek):
    Processing Chunk  5/15 ( 33%, est:    0 min:03 sek):
    Processing Chunk  6/15 ( 40%, est:    0 min:03 sek):
    Processing Chunk  7/15 ( 47%, est:    0 min:02 sek):
    Processing Chunk  8/15 ( 53%, est:    0 min:02 sek):
    Processing Chunk  9/15 ( 60%, est:    0 min:02 sek):
    Processing Chunk 10/15 ( 67%, est:    0 min:02 sek):
    Processing Chunk 11/15 ( 73%, est:    0 min:01 sek):
    Processing Chunk 12/15 ( 80%, est:    0 min:01 sek):
    Processing Chunk 13/15 ( 87%, est:    0 min:01 sek):
    Processing Chunk 14/15 ( 93%, est:    0 min:01 sek):
    Processing Chunk 15/15 (100%, est:    0 min:00 sek):

    Result for file: C_Easy1_noise02.mat
             |U2|Spks|  NO|  O| Det|FP|FN|FN NO|FN O|  TP|TP NO|TP O|Cl|Cl NO|Cl O|DetEO|DetE|EO| E |
    ---------|--------------------------------------------------------------------------------------|
    GT Unit 1| 1 1198 1067 131 1193  0  5     0    5 1193  1067  126  0     0    0     5    5  5  5 |
    GT Unit 2| 2 1128 1014 114 1123  0  5     0    5 1123  1014  109  0     0    0     5    5  5  5 |
    GT Unit 3| 3 1148 1031 117 1148  0  0     0    0 1148  1031  117  0     0    0     0    0  0  0 |
    ---------|--------------------------------------------------------------------------------------|
    Total    | 6 3474 3112 362 3464  0 10     0   10 3464  3112  352  0     0    0    10   10 10 10 |
	

Compare the results in the last row of the last table with the 4th row of the 
first table. The last table contains the evaluation of the BOTM sorter. Keep
in mind that we "cheated" by using the correct templates and correct noise
statistics. This in a way supervised BOTM produces 10 errors in total (detection 
and classification together). All of those 10 errors are not detected spikes in
an overlap.
In the original 2004 publication (the sorting here was unsupervised!) Quiroga et
al report 12 classification errors on the 2678 non overlapping spikes and 
714+275+10=999 detection errors on all spikes.
The first table contains a lot of "-1" entries. These numbers are not available
from the 2004 publication.

It is hard to compare the result from an unsupervised algorithm and this of an 
supervised one, that uses the ground truth for its initialisation.
However, if we repeat the analysis with the filter spike sorter (BOTM) 
and use the output of wave_clus instead of the real ground truth to calculate 
templates and noise statistics, we arrive at a similar result (around 50 
detection errors and 7 overlap errors). That means the BOTM is able to reduce the
error rate dramatically. This analysis is not included in this package.