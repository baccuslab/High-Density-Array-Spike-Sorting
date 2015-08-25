function [Q BIGTABLE] = E10_QuirogaPerformanceIn2004Paper()
    % Numbers taken from 
    % "Unsupervised spike detection and sorting with 
    %  wavelets and superparamagnetic clustering", Quiroga et al, 2004.
    %
    % The first table gives, "detection errors"
    % the second table the classification errors

    Q{1}.name = 'Easy 1';
    %                          NOISE   #SP  OVPs FN  FNO  FP
    Q{1}.DetectionErrors = [ 0.05 3514 (785) 17 (193) 711
                                0.10 3522 (769) 2 (177) 57
                                0.15 3477 (784) 145 (215) 14
                                0.20 3474 (796) 714 (275) 10];
    %                                |          SPC         |  K-Means
    %                       NOISE #Sp Wav PCA SpShap FeatSet   Wav PCA                         
    Q{1}.ClassificationErrors = [0.05 2729 1 1 0 863 0 0
                            0.10 2753 5 17 0 833 0 0
                            0.15 2693 5 19 0 2015 0 0
                            0.20 2678 12 130 24 614 17 17
                            0.25 2586 64 911 266 1265  69 68
                            0.30 2629 276 1913 838 1699  177 220
                            0.35 2702 483 1926  1424  1958  308 515
                            0.40 2645 741 1738 1738  1977  930 733];

    Q{2}.name = 'Easy 2';
    %                          NOISE   #SP  OVPs FN  FNO  FP
    Q{2}.DetectionErrors = [ 0.05 3410 (791) 0 (174) 0
                                0.10 3520 (826) 0 (191) 2
                                0.15 3411 (763) 10 (173) 1
                                0.20 3526 (811) 376 (256) 5];
    %                                |          SPC         |  K-Means
    %                       NOISE #Sp Wav PCA SpShap FeatSet   Wav PCA                          
    Q{2}.ClassificationErrors = [0.05 2619 3 4 2 502 0 0
                            0.10 2694 10 704 59 1893  2 53
                            0.15 2648 45 1732  1054  2199  31 336
                            0.20 2715 306 1791  2253 2199  154 740];

    Q{3}.name = 'Difficult 1';
    %                          NOISE   #SP  OVPs FN  FNO  FP
    Q{3}.DetectionErrors = [ 0.05 3383 (767) 1 (210) 63
                                0.10 3448 (810) 0 (191) 10
                                0.15 3472 (812) 8 (203) 6
                                0.20 3414 (790) 184 (219) 2];
    %                                |          SPC         |  K-Means
    %                       NOISE #Sp Wav PCA SpShap FeatSet   Wav PCA                          
    Q{3}.ClassificationErrors = [0.05 2616 0 7 3 619 0 1
                            0.10 2638 41 1781 794 1930  850 184
                            0.15 2660 81 1748  2131  2150  859 848
                            0.20 2624 651 1711  2449  2185  874 1170];

    Q{4}.name = 'Difficult 2';
    %                          NOISE   #SP  OVPs FN  FNO  FP
    Q{4}.DetectionErrors = [ 0.05 3364 (829) 0 (182) 1
                                0.10 3462 (720) 0 (152) 5
                                0.15 3440 (809) 3 (186) 4
                                0.20 3493 (777) 262 (228) 2];
    %                                |          SPC         |  K-Means
    %                       NOISE #Sp Wav PCA SpShap FeatSet   Wav PCA                          
    Q{4}.ClassificationErrors = [0.05 2535 1 1310 24 1809 686 212
                            0.10 2742 8 946  970  1987  271 579
                            0.15 2631 443 1716  1709  2259  546 746
                            0.20 2716 1462  1732  1732  1867  872 1004];


    for i=1:4
        nRows = size(Q{i}.ClassificationErrors,1);
        Q{i}.Table = -1*ones(nRows, 7);

        for j=1:nRows
            Q{i}.rowLabel{j} = sprintf('Noise %.2f', Q{i}.ClassificationErrors(j,1));        
            if j <= size(Q{i}.DetectionErrors,1)
                c = 1;
                Q{i}.colLabel{c} = '#Spks'; 
                Q{i}.Table(j,c) = Q{i}.DetectionErrors(j,2); c=c+1;
                Q{i}.colLabel{c} = '#Ovps';
                Q{i}.Table(j,c) = Q{i}.DetectionErrors(j,3); c=c+1;
                Q{i}.colLabel{c} = '#FN';
                Q{i}.Table(j,c) = Q{i}.DetectionErrors(j,4); c=c+1;
                Q{i}.colLabel{c} = '#FNOv';
                Q{i}.Table(j,c) = Q{i}.DetectionErrors(j,5); c=c+1;
                Q{i}.colLabel{c} = '#FP';
                Q{i}.Table(j,c) = Q{i}.DetectionErrors(j,6); c=c+1;
            end
            c=6;
            Q{i}.colLabel{c} = 'clas.Spks';
            Q{i}.Table(j,c) = Q{i}.ClassificationErrors(j,2); c=c+1;
            Q{i}.colLabel{c} = 'clas.Err';
            Q{i}.Table(j,c) = Q{i}.ClassificationErrors(j,3); c=c+1;        
        end
    end
    BIGTABLE = [];
    for i=1:4
        BIGTABLE = [BIGTABLE; Q{i}.Table];
    end    
    % for i=1:4
    %     mysort.plot.printTable(Q{i}.Table(:,1:7), 'colLabel', Q{i}.colLabel,...
    %             'rowLabel', Q{i}.rowLabel(:), 'topLeftLabel', Q{i}.name, 'printColSum', 1);
    % end        
