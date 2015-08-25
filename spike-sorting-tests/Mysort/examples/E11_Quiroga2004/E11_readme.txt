Here, we illustrate the systematic evaluation of the sorter on all data 
files in the Quiroga 2004 benchmark. The sorter is run with the default
values on all files. Then a table is printed, that contrasts the original
published error rates with the newly obtained:

				 |QDet|QCla|QTot|BDet|BCla|BTot |
	-------------------------|------------------------------|
	C_Easy1_noise005.mat     | 921    1  922|  14    1   15 |
	C_Easy1_noise01.mat      | 236    5  241|   3    0    3 |
	C_Easy1_noise015.mat     | 374    5  379|   7    0    7 |
	C_Easy1_noise02.mat      | 999   12 1011|  10    0   10 |
	-------------------------|------------------------------|
	C_Easy2_noise005.mat     | 174    3  177|   8    9   17 |
	C_Easy2_noise01.mat      | 193   10  203|   8    8   16 |
	C_Easy2_noise015.mat     | 184   45  229|  12   12   24 |
	C_Easy2_noise02.mat      | 637  306  943|  17   38   55 |
	-------------------------|------------------------------|
	C_Difficult1_noise005.mat| 274    0  274|   3   58   61 |
	C_Difficult1_noise01.mat | 201   41  242|  11  109  120 |
	C_Difficult1_noise015.mat| 217   81  298|  17  117  134 |
	C_Difficult1_noise02.mat | 405  651 1056|  13  212  225 |
	-------------------------|------------------------------|
	C_Difficult2_noise005.mat| 183    1  184|   9   20   29 |
	C_Difficult2_noise01.mat | 157    8  165|  19   26   45 |
	C_Difficult2_noise015.mat| 193  443  636|  18   69   87 |
	C_Difficult2_noise02.mat | 492 1462 1954|  30  233  263 |
	-------------------------|------------------------------|
	Total                    |5840 3074 8914| 199  912 1111 |

(Revision: 5260, TU-B)

Note the dramatic decrease of the total error number (8914 to 1570). But 
again keep in mind that we used the ground truth information to set up 
the BOTM sorter. The result using the output of wave_clus for the 
initialisation of the bss is similar though (not shown here).


Update 22.Mai 2012:
                             |QDet|QCla|QTot|BDet|BCla|BTot |
    -------------------------|------------------------------|
    C_Easy1_noise005.mat     | 921    1  922|   5    4    9 |
    C_Easy1_noise01.mat      | 236    5  241|   3    2    5 |
    C_Easy1_noise015.mat     | 374    5  379|   5    1    6 |
    C_Easy1_noise02.mat      | 999   12 1011|   4    3    7 |
    -------------------------|------------------------------|
    C_Easy2_noise005.mat     | 174    3  177|   3    0    3 |
    C_Easy2_noise01.mat      | 193   10  203|   5    0    5 |
    C_Easy2_noise015.mat     | 184   45  229|   5    3    8 |
    C_Easy2_noise02.mat      | 637  306  943|   7    2    9 |
    -------------------------|------------------------------|
    C_Difficult1_noise005.mat| 274    0  274|   7   17   24 |
    C_Difficult1_noise01.mat | 201   41  242|  17   14   31 |
    C_Difficult1_noise015.mat| 217   81  298|  10   18   28 |
    C_Difficult1_noise02.mat | 405  651 1056|  20   15   35 |
    -------------------------|------------------------------|
    C_Difficult2_noise005.mat| 183    1  184|   7   10   17 |
    C_Difficult2_noise01.mat | 157    8  165|   7   12   19 |
    C_Difficult2_noise015.mat| 193  443  636|   3   11   14 |
    C_Difficult2_noise02.mat | 492 1462 1954|   8   20   28 |
    -------------------------|------------------------------|
    Total                    |5840 3074 8914| 116  132  248 |

(Revision: 24049, ETH)