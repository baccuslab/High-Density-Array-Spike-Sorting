                         |QDet|QCla|QTot|BDet|BCla|BTot |
-------------------------|------------------------------|
C_Easy1_noise005.mat     | 921    1  922| 223    1  224 |
C_Easy1_noise01.mat      | 236    5  241| 202    2  204 |
C_Easy1_noise015.mat     | 374    5  379| 233    1  234 |
C_Easy1_noise02.mat      | 999   12 1011| 198    5  203 |
-------------------------|------------------------------|
C_Easy2_noise005.mat     | 174    3  177| 198   13  211 |
C_Easy2_noise01.mat      | 193   10  203| 216   14  230 |
C_Easy2_noise015.mat     | 184   45  229| 186   20  206 |
C_Easy2_noise02.mat      | 637  306  943| 310   86  396 |
-------------------------|------------------------------|
C_Difficult1_noise005.mat| 274    0  274| 209   12  221 |
C_Difficult1_noise01.mat | 201   41  242| 184   38  222 |
C_Difficult1_noise015.mat| 217   81  298| 197  113  310 |
C_Difficult1_noise02.mat | 405  651 1056| 176  259  435 |
-------------------------|------------------------------|
C_Difficult2_noise005.mat| 183    1  184| 201   24  225 |
C_Difficult2_noise01.mat | 157    8  165| 158   29  187 |
C_Difficult2_noise015.mat| 193  443  636| 204  138  342 |
C_Difficult2_noise02.mat | 492 1462 1954| 221  371  592 |
-------------------------|------------------------------|
Total                    |5840 3074 8914|3316 1126 4442 |

EpochLength 10
                         |QDet|QCla|QTot|BDet|BCla|BTot |
-------------------------|------------------------------|
C_Easy1_noise005.mat     | 921    1  922|  89    1   90 |
C_Easy1_noise01.mat      | 236    5  241|  59    2   61 |
C_Easy1_noise015.mat     | 374    5  379|  84    1   85 |
C_Easy1_noise02.mat      | 999   12 1011|  81    5   86 |
-------------------------|------------------------------|
C_Easy2_noise005.mat     | 174    3  177|  60   19   79 |
C_Easy2_noise01.mat      | 193   10  203|  50   21   71 |
C_Easy2_noise015.mat     | 184   45  229|  54   38   92 |
C_Easy2_noise02.mat      | 637  306  943| 167   93  260 |
-------------------------|------------------------------|
C_Difficult1_noise005.mat| 274    0  274|  44   17   61 |
C_Difficult1_noise01.mat | 201   41  242|  49   49   98 |
C_Difficult1_noise015.mat| 217   81  298|  55  127  182 |
C_Difficult1_noise02.mat | 405  651 1056|  54  275  329 |
-------------------------|------------------------------|
C_Difficult2_noise005.mat| 183    1  184|  53   30   83 |
C_Difficult2_noise01.mat | 157    8  165|  43   37   80 |
C_Difficult2_noise015.mat| 193  443  636|  65  148  213 |
C_Difficult2_noise02.mat | 492 1462 1954|  95  388  483 |
-------------------------|------------------------------|
Total                    |5840 3074 8914|1102 1251 2353 |


                         |QDet|QCla|QTot|BDet|BCla|BTot |
-------------------------|------------------------------|
C_Easy1_noise005.mat     | 921    1  922|  87    0   87 |
C_Easy1_noise01.mat      | 236    5  241|  68    0   68 |
C_Easy1_noise015.mat     | 374    5  379|  90    1   91 |
C_Easy1_noise02.mat      | 999   12 1011|  84    0   84 |
-------------------------|------------------------------|
C_Easy2_noise005.mat     | 174    3  177|  63    7   70 |
C_Easy2_noise01.mat      | 193   10  203|  50    9   59 |
C_Easy2_noise015.mat     | 184   45  229|  53   14   67 |
C_Easy2_noise02.mat      | 637  306  943|  68   30   98 |
-------------------------|------------------------------|
C_Difficult1_noise005.mat| 274    0  274|  58   93  151 |
C_Difficult1_noise01.mat | 201   41  242|  59  138  197 |
C_Difficult1_noise015.mat| 217   81  298|  59  126  185 |
C_Difficult1_noise02.mat | 405  651 1056|  63  224  287 |
-------------------------|------------------------------|
C_Difficult2_noise005.mat| 183    1  184|  52   25   77 |
C_Difficult2_noise01.mat | 157    8  165|  43   24   67 |
C_Difficult2_noise015.mat| 193  443  636|  65   60  125 |
C_Difficult2_noise02.mat | 492 1462 1954|  66  227  293 |
-------------------------|------------------------------|
Total                    |5840 3074 8914|1028  978 2006 |


    alpha = .5;
    Covest.CCol = mysort.noise.Cte2Ccol((1-alpha)*GT.C + alpha*diag(diag(GT.C)), 1);
%     NE = mysort.util.NoiseEstimator(GT.C, Tf);
%     botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 3, 'spikePrior', .01, 'chunk_size', 500000);
    botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 1, 'spikePrior', .0001, 'useSIC', false, 'minEpochLength', 5, 'chunk_size', 500000);

                         |QDet|QCla|QTot|BDet|BCla|BTot |
-------------------------|------------------------------|
C_Easy1_noise005.mat     | 921    1  922| 810    8  818 |
C_Easy1_noise01.mat      | 236    5  241| 710   14  724 |
C_Easy1_noise015.mat     | 374    5  379| 586   13  599 |
C_Easy1_noise02.mat      | 999   12 1011| 406   17  423 |
-------------------------|------------------------------|
C_Easy2_noise005.mat     | 174    3  177|  47   20   67 |
C_Easy2_noise01.mat      | 193   10  203|  35   23   58 |
C_Easy2_noise015.mat     | 184   45  229|  39   37   76 |
C_Easy2_noise02.mat      | 637  306  943|  82   94  176 |
-------------------------|------------------------------|
C_Difficult1_noise005.mat| 274    0  274|  26   21   47 |
C_Difficult1_noise01.mat | 201   41  242|  30   53   83 |
C_Difficult1_noise015.mat| 217   81  298|  37  135  172 |
C_Difficult1_noise02.mat | 405  651 1056|  36  281  317 |
-------------------------|------------------------------|
C_Difficult2_noise005.mat| 183    1  184|  39   36   75 |
C_Difficult2_noise01.mat | 157    8  165|  28   40   68 |
C_Difficult2_noise015.mat| 193  443  636|  42  155  197 |
C_Difficult2_noise02.mat | 492 1462 1954| 110  386  496 |
-------------------------|------------------------------|
Total                    |5840 3074 8914|3063 1333 4396 |


    alpha = .5;
    Covest.CCol = mysort.noise.Cte2Ccol((1-alpha)*GT.C + alpha*diag(diag(GT.C)), 1);
%     NE = mysort.util.NoiseEstimator(GT.C, Tf);
%     botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 3, 'spikePrior', .01, 'chunk_size', 500000);
    botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 1, 'spikePrior', .0001, 'useSIC', false, 'minEpochLength', 8, 'chunk_size', 500000);

                         |QDet|QCla|QTot|BDet|BCla|BTot |
-------------------------|------------------------------|
C_Easy1_noise005.mat     | 921    1  922|  76    2   78 |
C_Easy1_noise01.mat      | 236    5  241|  45    3   48 |
C_Easy1_noise015.mat     | 374    5  379|  59    5   64 |
C_Easy1_noise02.mat      | 999   12 1011|  61    9   70 |
-------------------------|------------------------------|
C_Easy2_noise005.mat     | 174    3  177|  51   19   70 |
C_Easy2_noise01.mat      | 193   10  203|  39   21   60 |
C_Easy2_noise015.mat     | 184   45  229|  43   36   79 |
C_Easy2_noise02.mat      | 637  306  943|  85   93  178 |
-------------------------|------------------------------|
C_Difficult1_noise005.mat| 274    0  274|  37   17   54 |
C_Difficult1_noise01.mat | 201   41  242|  39   50   89 |
C_Difficult1_noise015.mat| 217   81  298|  46  132  178 |
C_Difficult1_noise02.mat | 405  651 1056|  45  278  323 |
-------------------------|------------------------------|
C_Difficult2_noise005.mat| 183    1  184|  45   34   79 |
C_Difficult2_noise01.mat | 157    8  165|  33   38   71 |
C_Difficult2_noise015.mat| 193  443  636|  51  153  204 |
C_Difficult2_noise02.mat | 492 1462 1954| 104  386  490 |
-------------------------|------------------------------|
Total                    |5840 3074 8914| 859 1276 2135 |


    alpha = .5;
    Covest.CCol = mysort.noise.Cte2Ccol((1-alpha)*GT.C + alpha*diag(diag(GT.C)), 1);
%     NE = mysort.util.NoiseEstimator(GT.C, Tf);
%     botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 3, 'spikePrior', .01, 'chunk_size', 500000);
    botm = mysort.sorters.BOTM(Covest, Tf, GT.templates, 'upsample', 3, 'spikePrior', .0001, 'useSIC', true, 'minEpochLength', 8, 'chunk_size', 500000);
                         |QDet|QCla|QTot|BDet|BCla|BTot |
-------------------------|------------------------------|
C_Easy1_noise005.mat     | 921    1  922|  33    1   34 |
C_Easy1_noise01.mat      | 236    5  241|  17    3   20 |
C_Easy1_noise015.mat     | 374    5  379|  24    4   28 |
C_Easy1_noise02.mat      | 999   12 1011|  30    6   36 |
-------------------------|------------------------------|
C_Easy2_noise005.mat     | 174    3  177|   8   23   31 |
C_Easy2_noise01.mat      | 193   10  203|   9   22   31 |
C_Easy2_noise015.mat     | 184   45  229|  15   38   53 |
C_Easy2_noise02.mat      | 637  306  943|  62   86  148 |
-------------------------|------------------------------|
C_Difficult1_noise005.mat| 274    0  274|   2   15   17 |
C_Difficult1_noise01.mat | 201   41  242|  11   19   30 |
C_Difficult1_noise015.mat| 217   81  298|  15   67   82 |
C_Difficult1_noise02.mat | 405  651 1056|  22  194  216 |
-------------------------|------------------------------|
C_Difficult2_noise005.mat| 183    1  184|   9   36   45 |
C_Difficult2_noise01.mat | 157    8  165|  12   35   47 |
C_Difficult2_noise015.mat| 193  443  636|  23  141  164 |
C_Difficult2_noise02.mat | 492 1462 1954|  90  360  450 |
-------------------------|------------------------------|
Total                    |5840 3074 8914| 382 1050 1432 |
