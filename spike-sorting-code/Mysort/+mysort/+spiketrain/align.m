
function R = align(St1, St2, maxJitter, maxShift, maxOverlapDist, St1IDs, St2IDs)
% Alignment and assignment of spike trains used for evaluation purposes
% This function takes two sets of spike trains (point processes) and
% computes all pairwise comparisons between the two sets. Then, starting
% with the best matching pair, it assigns spike trains to each other in a
% greedy fashion.
% Once a spike train from on set is assigned to a spike
% train in the other set, it cannot be assigned to another spike train
% anymore.
% If St1 is the ground truth (e.g. from a spike sorting benchmark) and St2
% is the set of sorted spike trains, this function can give you an
% evaluation of how well the sorting was done. To visualize the resulting 
% assignments use:
%
%     R = mysort.spiketrain.alignGDFs(gdf1, gdf2, 10, 10, 10);
%     mysort.plot.printEvaluationTable(R, 'ignoreEmptySorted', 0);
%     mysort.plot.alignment(R)
% 
% Inputs:
%   St1 - cell array containing spike trains. Each spike train is an array
%         of spike times. St1 = {[1 10 100], [50 150 1000]};
%         In this function St1 is thought to be the ground truth but in
%         principle St1 and St2 can be exchanged.
%   St2 - second set of spike trains, same format as St1.
%  maxJitter - maximum time difference between two spikes in St1 and St2 so
%              that they are still assumed to be the same spike
%  maxShift  - maximum time shift of a whole spike train of St2. During the
%              matching process, whole spike trains of St2 are shifted to
%              see if there is a consisten time lag between spike trains in
%              St2 and St1. This could arise by the time in the ground
%              truth (St1) being defined as the beginning of a spike
%              waveform, while the sorted spike trains (St2) give the time
%              of the peak of the waveforms.
%  maxOverlapDist - Spiketrains in St1 are compared to estimate whether
%                   there are overlapping spikes. maxOverlapDist is the
%                   maximum time difference between to spikes in to
%                   different spike trains of St1 so that they are called
%                   overlaps (remember, St1 is the ground truth, so it
%                   defines what is an overlap and what not).
%
% Outputs:
%    R  - big structure with information about the assignement process and
%         the final results
%
    R.St1IDs = 1:length(St1);
    if exist('St1IDs', 'var')
        R.St1IDs = St1IDs(:)';
    end
    R.St2IDs = 1:length(St2);
    if exist('St2IDs', 'var')
        R.St2IDs = St2IDs(:)';
    end    
    R.nSt1 = length(St1);
    R.nSt2 = length(St2);
    def = mysort.util.defs();
    % Make sure the spike times in the spike trains are sorted
    for s1=1:R.nSt1;  St1{s1} = sort(St1{s1}); end
    for s2=1:R.nSt2;  St2{s2} = sort(St2{s2}); end
    
    % Compute all pairwise spike train comparisons between a spike train
    % from the ground truth (st1) and a spike train from st2 (sorted 
    % spike trains). 
    R.RsingleComp = mysort.spiketrain.compare(St1, St2, maxJitter, maxShift);
    
    % For every ground truth spike train, find a sorted spike train, that
    % fits best. This is the one, that has the most true positives (TP) 
    % with this ground truth spike train.
    R = computeSpikeTrainAssignment(R, St1, St2);
    
    % Recalculate the comparison and the assignments with the shifted spike
    % trains
    if maxShift > 0
        % Shift the spike trains according to the alignment
        St2 = shiftSpikeTrains(R, St2);
        R = mysort.spiketrain.align(St1, St2, maxJitter, 0, maxOverlapDist, R.St1IDs, R.St2IDs);
        return
    end
    
    % Check every ground truth spike, whether there is another one nearby
    % as an overlap.
    [R.O R.nO] = mysort.spiketrain.checkForOverlaps(St1, maxOverlapDist);
    
    % Label every spike that can be associated with another spike
    R = computeLabels(R, St1, St2, maxJitter);
    
    % Now we can label everything that has no label. These are the FPs and 
    % FNs respectively
    R = labelFP_FN(R);
    
    % Compute the statistics
    R = computeStats(R);
    R.RsingleComp.readme = 'This struct contains the results of the purely pairwise comparison between the ground truth and sorted spike trains.';

    R = buildTable(R);
    R.St1 = St1;
    R.St2 = St2;
    
    R = buildSimilarityMatrix(R);
    R.readme = 'This struct contains the information for every spike of the ground truth spike trains and the sorted spike trains, whether it is a true positive, an error (FN, FP) and whether it participates in an overlap. For this, pairwise assignements are were computed and stored. Furthermore, if a spike from a ground truth spike train and a sorted spike train were matched, then this alignment is also stored.';
    


    %----------------------------------------------------------------------
    function e = sortedErr2associatedGTErr(errs)
        % This function reorders the errors associated with sorted spike
        % trains in a way that they are ordered according to the associated
        % ground truth spike trains
        e = errs(R.k2f(R.k2f>0));
        if R.nSt1 <= R.nSt2
            return
        else
            tmp = zeros(1, R.nSt1);
            tmp(R.k2f>0) = errs(R.k2f(R.k2f>0));
            e = tmp;
        end
    end
 
    %----------------------------------------------------------------------
    function R = labelFP_FN(R)
        % Spikes of the Ground truth that have no label are false negatives
        for i=1:R.nSt1
            R.spikeLabel1{i}(R.spikeLabel1{i}==-1 & R.O{i}) = def.fno;
            R.spikeLabel1{i}(R.spikeLabel1{i}==-1) = def.fn;
        end
        % Spikes of the sorted spike trains that have no label are false
        % positives
        for j=1:R.nSt2
            R.spikeLabel2{j}(R.spikeLabel2{j}==-1) = def.fp;
        end
    end
    
    %----------------------------------------------------------------------
    function St2 = shiftSpikeTrains(R, St2)
        for i=1:R.nSt1
            if R.k2f(i)>0
                j = R.k2f(i);
                St2{j} = St2{j} + R.RsingleComp.shift(i,j); 
            end
        end
    end

    %----------------------------------------------------------------------
    function R = computeSpikeTrainAssignment(R, St1, St2)
        %nAssignments = min(R.nSt1, R.nSt2);
        TP = R.RsingleComp.TP;
        blockedRows = [];
        blockedCols = [];
        R.k2f = ones(1, R.nSt1)*-1;
        R.f2k = ones(1, R.nSt2)*-1;
        R.spikeLabel1 = {};
        R.spikeLabel2 = {};        
        for i=1:length(St1)
            R.spikeLabel1{i} = ones(1,length(St1{i}))*-1;
        end
        for j=1:length(St2)
            R.spikeLabel2{j} = ones(1,length(St2{j}))*-1;
        end       
        count = 0;
        while count < R.nSt1*R.nSt2
            [i,j,v] = mysort.util.matrixArgMax(TP);

            if ~any(i==blockedRows) && ~any(j==blockedCols)
                % This means we are in one of the "correct" assignments
                % between two spike trains that fit best together
                blockedRows = [blockedRows i];
                blockedCols = [blockedCols j];
                R.k2f(i) = j;
                R.f2k(j) = i;          
            end
            if length(blockedRows) >= R.nSt1 || length(blockedCols)>=R.nSt2
                return
            end
            TP(i,j) = -1;
            count = count+1;
        end
    end

    %----------------------------------------------------------------------
    function R = computeLabels(R, St1, St2, maxJitter)
        % Compute the label for the "correct" associations first
        for i=1:R.nSt1
            j = R.k2f(i);
            if j>0
                associationFlag = 1;
                [R.ALI{i,j} R.spikeLabel1{i} R.spikeLabel2{j}] = ...
                    computeAlignment(St1{i}, St2{j}, associationFlag, ...
                        R.spikeLabel1{i}, R.spikeLabel2{j}, maxJitter, R.O{i});
            end
        end
        % Now compute all that are "wrongly" associated
        TP = R.RsingleComp.TP;
        
        count = 0;
        while count < R.nSt1*R.nSt2
            count = count+1;
            [i,j,v] = mysort.util.matrixArgMax(TP);
            TP(i,j) = -1;
            if R.k2f(i)>0 && R.f2k(j) == i
                % These two spike trains are assigned and already computed!
                
                continue
            elseif R.k2f(i)>0 && R.f2k(j)>0
                % This is an assignment between two spike trains that do 
                % not really belong together. BOTH spike trains are
                % already in an association with another spike train (but
                % not with each other!)

            elseif R.k2f(i)>0
                % This is an assignment between two spike trains that do 
                % not really belong together. 
                % The GROUNDTRUTH spike train is already blocked.

            elseif R.f2k(j)>0
                % This is an assignment between two spike trains that do 
                % not really belong together. 
                % The SORTED spike train is already blocked.

            else
                % This means that neither of the spike trains is associated
                % should not happen!
                error('Why?!');         
            end
            associationFlag = 0;
            [R.ALI{i,j} R.spikeLabel1{i} R.spikeLabel2{j}] = ...
                computeAlignment(St1{i}, St2{j}, associationFlag, ...
                    R.spikeLabel1{i}, R.spikeLabel2{j}, maxJitter, R.O{i});
        end              
    end

    %----------------------------------------------------------------------
    function [ali, label1, label2] = computeAlignment(st1, st2, ...
            associationFlag, label1, label2, maxJitter, overlap)
        % Compute for this pair of spike trains the one to one
        % relationships between their spikes. If a spike has already a
        % label, ignore it.
        idxI = 1;
        idxJ = 1;
        ali = [];
        while idxI <= length(st1) && idxJ <= length(st2)
            if label1(idxI) ~= -1
                idxI = idxI +1;
                continue
            end
            
            if label2(idxJ) ~= -1
                idxJ = idxJ +1;
                continue
            end
            
            if st1(idxI) <= st2(idxJ) + maxJitter && ...
               st1(idxI) >= st2(idxJ) - maxJitter
                ali = [ali; [idxI idxJ]];
                if associationFlag && ~overlap(idxI)
                    label1(idxI) = def.tp;
                    label2(idxJ) = def.tp;
                elseif associationFlag && overlap(idxI)
                    label1(idxI) = def.tpo;
                    label2(idxJ) = def.tpo;
                elseif ~associationFlag && ~overlap(idxI)
                    label1(idxI) = def.cl;
                    label2(idxJ) = def.cl;
                elseif ~associationFlag && overlap(idxI)
                    label1(idxI) = def.clo;
                    label2(idxJ) = def.clo;  
                else
                    error('Whoopsi');
                end
                idxI = idxI +1;
                idxJ = idxJ +1;
            elseif st1(idxI) < st2(idxJ)
                idxI = idxI +1;
            else
                idxJ = idxJ +1;
            end
        end
    end

    %----------------------------------------------------------------------
    function R = computeStats(R)
        R.nSP1   = zeros(1,R.nSt1);
        R.nSP1NO = zeros(1,R.nSt1);
        R.nSP1O  = zeros(1,R.nSt1);

        R.nTPNO = zeros(1,R.nSt1);
        R.nTPO  = zeros(1,R.nSt1);
        R.nFNNO = zeros(1,R.nSt1);
        R.nFNO  = zeros(1,R.nSt1);
        R.nCLNO = zeros(1,R.nSt1);
        R.nCLO  = zeros(1,R.nSt1);
        
        R.nSP2   = zeros(1,R.nSt2);
        R.nFP2   = zeros(1,R.nSt2);
        R.nCLNO2 = zeros(1,R.nSt2);
        R.nCLO2  = zeros(1,R.nSt2);
	
        for i=1:R.nSt1
            R.nSP1(i) = length(R.O{i});
            R.nSP1NO(i)= sum(R.O{i}==0);
            R.nSP1O(i)= sum(R.O{i});
            R.nTPNO(i) = sum(R.spikeLabel1{i}==def.tp);
            R.nTPO (i) = sum(R.spikeLabel1{i}==def.tpo);
            R.nFNNO(i) = sum(R.spikeLabel1{i}==def.fn);
            R.nFNO (i) = sum(R.spikeLabel1{i}==def.fno);
            R.nCLNO(i) = sum(R.spikeLabel1{i}==def.cl);
            R.nCLO (i) = sum(R.spikeLabel1{i}==def.clo);
        end
        for j=1:R.nSt2
            R.nSP2(j)  = length(R.spikeLabel2{j});
            R.nFP2(j)   = sum(R.spikeLabel2{j}==def.fp);
            R.nCLNO2(j) = sum(R.spikeLabel2{j}==def.cl);
            R.nCLO2 (j) = sum(R.spikeLabel2{j}==def.clo);            
        end  
        
        R.nDetected = sortedErr2associatedGTErr(R.nSP2);
        R.nFP       = sortedErr2associatedGTErr(R.nFP2);
        
        
        R.nTP = R.nTPNO + R.nTPO;
        R.nFN = R.nFNNO + R.nFNO;
        R.nCL = R.nCLNO + R.nCLO;
        R.nCL2= R.nCLNO2 + R.nCLO2;
        R.nFA = sortedErr2associatedGTErr(R.nCL2);
    
        R.detErr   = sortedErr2associatedGTErr(R.nFP2) + R.nFNNO + R.nFNO;
        R.detErrNO = sortedErr2associatedGTErr(R.nFP2) + R.nFNNO;
        R.detErrO  = R.nFNO;
        
        R.totErr   = R.nCL + R.detErr;
        R.totErrNO = R.nCLNO + R.detErrNO;
        R.totErrO  = R.nCLO + R.detErrO;
        
    end
    %----------------------------------------------------------------------
    function R = buildSimilarityMatrix(R)
        N = repmat(R.nSP1(:), 1, R.nSt2);
        S = zeros(size(R.ALI));
        for i=1:size(R.ALI,1)
            for j=1:size(R.ALI,2);
                S(i,j) = size(R.ALI{i,j},1);
            end
        end
        R.similaritySt1St2 = S;
        R.similaritySt1St2_normalized = S./N;
        
        S = zeros(size(R.RsingleComp.ALI));
        for i=1:size(R.RsingleComp.ALI,1)
            for j=1:size(R.RsingleComp.ALI,2);
                S(i,j) = size(R.RsingleComp.ALI{i,j},1);
            end
        end
        R.similarityBeforeAssignmentSt1St2 = S;     
        R.similarityBeforeAssignmentSt1St2_normalized = S./N;  
    end
    
    %----------------------------------------------------------------------
    function R = buildTable(R)
        assignedIDs = R.k2f(1:R.nSt1);
        idx = assignedIDs~=-1;
        assignedIDs(idx) = R.St2IDs(assignedIDs(idx));
        R.table = [R.St1IDs
                   assignedIDs
                   R.nSP1
                   R.nSP1NO
                   R.nSP1O
                   R.nDetected
                   R.nFA
                   R.nFP
                   R.nFN
                   R.nFNNO
                   R.nFNO
                   R.nTP
                   R.nTPNO
                   R.nTPO
                   R.nCL
                   R.nCLNO
                   R.nCLO
                   R.detErrO
                   R.detErr
                   R.totErrO
                   R.totErr]';
        R.tableLabelLong = {'True Unit', 'assigned Unit', ...
                        'True Spikes', 'Non-Overlaps', 'Overlaps', 'Detected',...
                        'Falsey Assigned from other Neuron', 'False Positive', 'False Negative', 'False Negative Non-Overlaps', ...
                        'False Negative Overlap',...
                        'True Positive', 'True Positive Non-Overlaps', 'True Positive Overlaps',...
                        'Classification Errors', 'Classification Errors Non-Overlaps', 'Classification Errors Overlaps',...
                        'Detection Errors Overlaps', 'total Detection Errors', ...
                        'total Errors Overlaps', 'total Errors'};
        R.tableLabel = {'True Unit', 'ass.Unit', ...
                        'True Spikes', 'Non-Overlaps', 'Overlaps', 'Detected',...
                        'False Assigned', 'False Positive', 'False Negative', 'False Negative NO', ...
                        'False Negative O',...
                        'True Positive', 'True Positive NO', 'True Positive O',...
                        'Classification Errors', 'Classification Errors NO', 'Classification Errors O',...
                        'Det Err O','tot Det Err',  ...
                        'tot Err O', 'tot Err'}; 
        R.tableLabelShort = {'U1', 'U2',...
                            'Spks', 'NO', 'O', 'Det', 'FA',...
                            'FP', 'FN', 'FN NO', ...
                            'FN O',...
                            'TP', 'TP NO', 'TP O',...
                            'Cl', 'Cl NO', 'Cl O',...
                            'DetEO', 'DetE', ...
                            'EO', 'E'}; 
        R.tableRowLabel = {};
        for i=1:R.nSt1
            R.tableRowLabel{i} = sprintf('GT Unit %d', R.table(i,1));
        end
        % Add sorted units that were not assigned to any GT unit
        notAssigned = find(R.f2k == -1);
        z = zeros(1, length(notAssigned));
        R.table = [R.table; [-1*ones(1, length(notAssigned))
                   R.St2IDs(notAssigned)
                   z
                   z
                   z
                   R.nSP2(notAssigned)
                   R.nCL2(notAssigned) % R.nFA
                   R.nSP2(notAssigned) - R.nCL2(notAssigned) %R.nFP
                   z
                   z
                   z
                   z
                   z
                   z
                   z % R.nCL2(notAssigned) % R.nCL
                   R.nCLNO2(notAssigned)
                   R.nCLO2(notAssigned)
                   z % R.detErrO
                   R.nSP2(notAssigned) - R.nCL2(notAssigned) %R.detErr
                   z %R.totErrO
                   z ]'];        %R.totErr
        for i=R.nSt1+1:size(R.table,1)
            R.tableRowLabel{i} = ' - - - ';
        end        
        checkTable();
    end

    function checkTable()
        % performs checks on the table if the numbers are consistent
        %assert(all(R.totErr == R.detErr + );
        assert(all(R.nDetected == R.nFA + R.nFP + R.nTP), 'Detected, FP and TP dont add up!');
        assert(all(R.nTP == R.nTPNO + R.nTPO), 'TP, TPNO and TPO dont add up!');
        assert(all(R.totErr - (R.nSP1 - R.nTP) - R.nFP == 0), 'total errors, true spikes, true positives and FP dont add up!');
        %assert(sum(R.nCL) == sum(R.nFA), 'Sum of classifcation errors must be the same as sum of false assignments!');
    end
end