els1 = hidens_get_all_electrodes(1);
els2 = hidens_get_all_electrodes(2);
date_ = date();
readme = 'This file contains the information returned by hidens_get_all_electrodes(x) to be available also under windows';
save('hidens_get_all_electrodes.h5', 'els1', 'els2', 'date_', 'readme', '-v7.3');