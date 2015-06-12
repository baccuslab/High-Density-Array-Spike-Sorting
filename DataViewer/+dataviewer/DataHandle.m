classdef DataHandle < mysort.util.DebuggableClass
    properties
        DBH
        GDFBuffer
        EGDFBuffer
    end
    properties (Constant = true)
        samplesPerSecond = 32000;
    end
    
    methods
        %------------------------------------------------------------------
        function self = DataHandle(varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            self.init(varargin{:});
        end

        %------------------------------------------------------------------
        function init(self, varargin)
            self.P.log_function('Trying to establish connection...');            
            self.DBH = db.munk.Postgres(varargin{:});
            self.P.log_function('... done.');
            self.GDFBuffer  = mysort.util.HashBufferedFunction(@(x) self.getGDF_(x));
            self.EGDFBuffer = mysort.util.HashBufferedFunction(@(x) self.getEGDF_(x));
        end
        %------------------------------------------------------------------
        function close(self)
            self.P.log_function('Closing connection...');            
            self.DBH.close();
            self.P.log_function('... done.');            
        end
        
        %------------------------------------------------------------------
        function execute(self, varargin)
            self.P.log_function('Executing SQL...');
            t1 = tic;
            self.DBH.execute(varargin{:});
            t = toc(t1);
            self.P.log_function('... done.');            
        end        
        %------------------------------------------------------------------
        function [R colNames colTypes] = query(self, varargin)
            self.P.log_function('Querying DB...');
            t1 = tic;
            [R colNames colTypes] = self.DBH.query(varargin{:});
            t = toc(t1);
            self.P.log_function(sprintf('... Query done (%.2fs).', t));
        end
        %------------------------------------------------------------------
        function R = getSrate(self)
            R = 1/self.samplesPerSecond;
        end
        %------------------------------------------------------------------
        function R = getSamplesPerSecond(self)
            R = self.samplesPerSecond;
        end
        
        %------------------------------------------------------------------
        function R = getMonkeys(self)
            qstr = 'SELECT e.subject_name FROM experiment AS e GROUP BY e.subject_name';            
            R = self.query(qstr);  
        end
        %------------------------------------------------------------------
        function R = getExperiments(self, varargin)
            P = self.perpareInputs([], varargin{:});
            if isempty(P.monkeys)
                if isempty(P.experiments)
                    qstr = 'SELECT experiment_id, name FROM experiment';
                else
                    qstr = 'SELECT experiment_id, name FROM experiment WHERE name=''%s''';
                    qstr = sprintf(qstr, P.experiments);
                end
            else
                if isempty(P.experiments)
                    qstr = 'SELECT experiment_id, name FROM experiment AS e WHERE subject_name=''%s''';
                    qstr = sprintf(qstr, P.monkeys);
                else
                    qstr = 'SELECT experiment_id, name FROM experiment AS e WHERE subject_name=''%s'' AND  name=''%s''';
                    qstr = sprintf(qstr, P.monkeys, P.experiments);
                end
            end
            R = self.query(qstr); 
        end
        %------------------------------------------------------------------
        function R = getTetrodes(self, varargin)
            P = self.perpareInputs([], varargin{:});
            tetrodeNrString = [];
            if ~isempty(P.tetrodeNRs)
                tetrodeNrString = [' AND te.name IN (' util.dlmstring(P.tetrodeNRs, ', ', 'T%02d') ') '];
            end
            suffix = ' AND SUBSTRING(te.name,1,1) IN (''T'') ORDER BY te.name';
            if ~isempty(P.experiments)
                qstr = ['SELECT te.channel_group_id, te.name FROM channel_group as te JOIN experiment as e ON (e.experiment_id = te.experiment_id) WHERE e.name = ''%s'' ' tetrodeNrString suffix];
                qstr = sprintf(qstr, P.experiments);
            elseif ~isempty(P.experimentIDs)
                qstr = ['SELECT te.channel_group_id, te.name FROM channel_group as te WHERE te.experiment_id=%d ' tetrodeNrString suffix];
                qstr = sprintf(qstr, P.experimentIDs);
            else
                error('not implemented');
            end
            R = self.query(qstr); 
        end
        %------------------------------------------------------------------
        function R = getChannels(self, varargin)
            P = self.perpareInputs([], varargin{:});
            if ~isempty(P.tetrodeIDs)
                qstr = ['SELECT c.channel_id, c.index FROM channel as c JOIN channel_group as te ON (te.channel_group_id = c.channel_group_id) WHERE te.channel_group_id IN (' util.dlmstring(P.tetrodeIDs) ') ORDER BY c.index'];
            else
                error('not implemented');
            end
            R = self.query(qstr, {'Double', 'Double'});
        end
        %------------------------------------------------------------------
        function R = getBlocks(self, varargin)
            P = self.perpareInputs([], varargin{:});
            select_str = 'SELECT b.block_id, b.name FROM block AS b JOIN experiment as e ON (e.experiment_id = b.experiment_id)';
            assert(length(P.blocks) < 2, 'Cannot select more than one block!');
            blockstr = [' b.name = ''' P.blocks ''''];
            if isempty(P.monkeys) && isempty(P.experiments) && isempty(P.experimentIDs) 
                if isempty(P.blocks)
                    qstr = 'SELECT block_id, name FROM block';
                else
                    qstr = ['SELECT block_id, name FROM block WHERE name = ''' P.blocks ''''];
                end
            elseif ~isempty(P.monkeys)
                if isempty(P.blocks)
                    qstr = [select_str ' WHERE e.subject_name=''%s'''];
                else
                    qstr = [select_str ' WHERE e.subject_name=''%s'' AND ' blockstr];
                end
                qstr = sprintf(qstr, P.monkeys);
            elseif ~isempty(P.experiments)
                if isempty(P.blocks)
                    qstr = [select_str ' WHERE e.name=''%s'''];
                else
                    qstr = [select_str ' WHERE e.name=''%s'' AND ' blockstr];
                end                
                qstr = sprintf(qstr, P.experiments);
            elseif ~isempty(P.experimentIDs)
                if isempty(P.blocks)
                    qstr = [select_str ' WHERE e.experiment_id=%d'];
                else
                    qstr = [select_str ' WHERE e.experiment_id=%d AND ' blockstr];
                end               

                qstr = sprintf(qstr, P.experimentIDs);
            else
                error('not implemented');
            end
            R = self.query(qstr); 
        end
        %------------------------------------------------------------------
        function R = getTrials(self, varargin)
            P.ok = [];
            P.rewarded = [];
            P.ignoreFileErrorTrials = 0;
            P.load = [1:4];
            P.ignoreLoadErrorTrials = 0;
            P = self.perpareInputs(P, varargin{:});
            
            % HACK TO GET THIS WORKING WITH NEW VERSION
            P.ok =[];
            P.ignoreFileErrorTrials = false;
            P.load = [];
            P.ignoreLoadErrorTrials = true;
            P.rewarded = [];
            
            if ~isempty(P.ok);       ok_str  = [' AND t.ok='  num2str(P.ok)       ' ']; else ok_str = ''; end
            if P.ignoreFileErrorTrials; ignore_str = ' AND t.read_error=0 '; else ignore_str=''; end
%             if ~isempty(P.load)
%                 load_str = ['AND (t.load IN (' util.dlmstring(P.load) ') '];
%                 if ~P.ignoreLoadErrorTrials; load_str = [load_str ' OR t.load IS NULL)']; else load_str = [load_str ') ']; end
%             else
%                 if ~P.ignoreLoadErrorTrials;
%                     load_str = 'AND t.load IS NULL ';
%                 else
%                     R = [];
%                     return
%                 end
%             end
            load_str = '';
            
            if ~isempty(P.rewarded)
                if P.rewarded;
                    rew_str = [' AND t.reward_time > 0 '];
                else
                    rew_str = [' AND t.reward_time IS NULL '];
                end
            else
                rew_str= '';
            end
            select_str = 'SELECT t.segment_id, t.index FROM segment AS t JOIN block as b ON (t.block_id = b.block_id)';
            if (~isempty(P.experiments) && ~isempty(P.blocks))
                qstr = [select_str ' JOIN experiment as e ON (e.experiment_id = b.experiment_id) WHERE e.name=''%s'' AND b.name=''%s'' ' ok_str rew_str ignore_str load_str ' ORDER BY t.index '];
                qstr = sprintf(qstr, P.experiments, P.blocks);
            elseif ~isempty(P.blockIDs)
                qstr = [select_str ' WHERE b.block_id IN (' util.dlmstring(P.blockIDs) ') ' ok_str rew_str ignore_str load_str ' ORDER BY t.index '];
            elseif ~isempty(P.experimentIDs)
                qstr = ['SELECT t.segment_id, t.index FROM segment as t JOIN block as b ON (b.block_id = t.block_id) JOIN experiment as e ON (e.experiment_id = b.experiment_id) WHERE e.experiment_id IN (' util.dlmstring(P.experimentIDs) ') ' ok_str rew_str ignore_str load_str ' ORDER BY t.index'];
            else
                error('not implemented');
            end
            R = self.query(qstr); 
        end
        
        %------------------------------------------------------------------
        function C = getCovarianceMatrix(self, varargin)
            P.kind = [];
            P = self.perpareInputs(P, varargin{:});
            if ~isempty(P.trialIDs)
                trialstring = ['and t.segment_id IN (' util.dlmstring(P.trialIDs) ')'];
                trial_join_string = ['JOIN segment as t on t.segment_id = cov.segment_id'];
            else
                trialstring = 'and cov.segment_id IS NULL';
                trial_join_string = '';
            end
            
            if isempty(P.kind)
                covkindstr = '';
            elseif strcmp(P.kind, 'NULL');
                covkindstr = 'and cov.kind_code IS NULL';
            else
                covkindstr = ['and cov.kind_code =''' P.kind ''''];
            end
            
            if length(P.analysisIDs) == 1
                qstr = ['SELECT c1.index, c2.index, xc.data FROM xcorr as xc JOIN covariance as cov ON cov.covariance_id = xc.covariance_id ' trial_join_string ' join analysis as a on a.analysis_id = cov.analysis_id join channel as c1 on c1.channel_id = xc.chan1_id join channel as c2 on c2.channel_id = xc.chan2_id where a.analysis_id = ' num2str(P.analysisIDs) ' ' trialstring ' ' covkindstr ' and c1.index <= c2.index order by c1.index, c2.index ASC'];
            else
                error('not implemented!')
            end
            R = self.query(qstr); 
            C = {};
            for i=1:size(R,1)
                if ~isempty(R{i,3})
                    C{R{i,1}+1, R{i,2}+1} = R{i,3};
                else
                    C{R{i,1}+1, R{i,2}+1} = [];
                end
            end
        end
        
        %------------------------------------------------------------------
        function xcovs = getCovarianceMatrixCellStruct(self, varargin)
            C_ = self.getCovarianceMatrix(varargin{:});
            xcovs = {};
            if isempty(C_); return; end
                
            Tf = (length(C_{1,1})-1)/2;
            for i = 1:size(C_,1)
                for j = 1:size(C_,2);
                    xcovs{i,j}.Tf = Tf;
                    xcovs{i,j}.xcov = C_{i,j};
                end
            end
        end
        
        %------------------------------------------------------------------
        function M = getSpikeTrain(self, varargin)
            P.alignOnEvent = [];
            P = self.perpareInputs(P, varargin{:});
            constraint_str = [];
            assert(~isempty(P.trialIDs), 'No trialID provided!');
            assert(~isempty(P.unitIDs), 'No unitID provided!');
            if ~isempty(P.from); constraint_str = [constraint_str ' AND s.sample>=' num2str(P.from)]; end
            if ~isempty(P.to);   constraint_str = [constraint_str ' AND s.sample<=' num2str(P.to)]; end 
            if isempty(P.alignOnEvent)
                qstr = ['SELECT s.segment_id, t.index, s.unit_id, u.index, s.sample FROM spike as s JOIN unit as u ON u.unit_id=s.unit_id JOIN segment as t ON t.segment_id = s.segment_id WHERE s.segment_id IN (' util.dlmstring(P.trialIDs) ') AND s.unit_id IN (' util.dlmstring(P.unitIDs) ') ' constraint_str ' ORDER BY s.segment_id, s.unit_id, s.sample ASC'];
            else
                assert(length(P.alignOnEvent)==1, 'Only on one event can be align!');
                qstr = ['SELECT s.segment_id, t.index, s.unit_id, u.index, s.sample-es.sample FROM spike as s JOIN unit as u ON u.unit_id=s.unit_id JOIN segment as t ON t.segment_id = s.segment_id JOIN eventsample as es ON es.segment_id = t.segment_id WHERE s.segment_id IN (' util.dlmstring(P.trialIDs) ') AND s.unit_id IN (' util.dlmstring(P.unitIDs) ') ' constraint_str ' AND es.event = ' num2str(P.alignOnEvent) ' ORDER BY s.segment_id, s.unit_id, s.sample ASC'];
            end
            R = self.query(qstr);
            M = cell2mat(R);
        end
        
        %------------------------------------------------------------------
        function R = getTrialLength(self, varargin)
            P = self.perpareInputs([], varargin{:});

            qstr = ['SELECT t.segment_id, sa.value FROM segment as t JOIN segment_attribute as sa ON (t.segment_id=sa.segment_id) WHERE sa.name=''length'' AND t.segment_id IN (' util.dlmstring(P.trialIDs) ') '...
                    'AND sa.value IS NOT NULL ORDER BY t.segment_id ASC']; 
            R = self.query(qstr);
            R = [cell2mat(R(:,1)) cellfun(@(x) str2double(x), R(:,2))];
        end
        
        %------------------------------------------------------------------
        function M = getEventSampleAndTrialLength(self, varargin)
            P = self.perpareInputs([], varargin{:});
            
            assert(~isempty(P.trialIDs), 'trialIDs need to be specified!');
            assert(~isempty(P.eventIDs), 'eventIDs need to be specified!');            
            qstr = ['Select t.segment_id, t.index, t.length, es.sample FROM segment as t JOIN eventsample as es ON (es.segment_id = t.segment_id) WHERE es.event IN (' util.dlmstring(P.eventIDs) ') AND t.segment_id IN (' util.dlmstring(P.trialIDs) ') ORDER BY t.index, es.sample ASC']; 
            M = cell2mat(self.query(qstr));
        end
        
        %------------------------------------------------------------------
        function [GDF null_trials eventSamples] = getGDF(self, varargin)
            P.alignOnEvent = false;
            P.cellForMultipleTrials = false;
            P = self.perpareInputs(P, varargin{:});
            % If multiple trials are selected the sql statement down there
            % gets insanely slow. fix this or live with this workaround:
            null_trials = [];
            spiketrain_container = cell(length(P.trialIDs),length(P.unitIDs));
            import matlabfilecentral.progressbar.progressbar
            progressbar(0, 0)         % Initialize/reset two bars
            progressbar('Trials', 'Units');
            P_ = P;
            for trial = 1:length(P.trialIDs)
                for unit = 1:length(P.unitIDs)
                    P_.segment_idIDs = P.trialIDs(trial);
                    P_.unitIDs  = P.unitIDs(unit);
                    
                    pHash = self.P2Hash(P_, 'unitIDs', 'trialIDs', 'from', 'to');            
                    GDF = self.GDFBuffer.call(pHash, P_);                    
                    spiketrain_container{trial, unit} = GDF;
                    progressbar([], unit/length(P.unitIDs));
                end
                progressbar((trial)/length(P.trialIDs), []);                
            end                
%             GDF = self.getGDF_(P);

            if P.alignOnEvent
                eventSamples = self.getEvents('trialIDs', P.trialIDs, 'eventIDs', P.eventIDs);
                eventSamples = cell2mat(eventSamples(:,[1:3 5]));
                if isempty(eventSamples)
                    null_trials = 1:length(P.trialIDs);
                else
                    for trial = 1:length(P.trialIDs)
                        idx = find(eventSamples(:,4) == P.trialIDs(trial), 1);
                        if ~isempty(idx) && eventSamples(idx,1)>0
                            for unit = 1:length(P.unitIDs)
                                if ~isempty(spiketrain_container{trial, unit})
                                    spiketrain_container{trial, unit}(:,2) = spiketrain_container{trial, unit}(:,2) - eventSamples(idx,2);
                                end
                            end
                        else
                            null_trials = [null_trials trial];
                            for unit = 1:length(P.unitIDs)
                                spiketrain_container{trial, unit} = [];
                            end
                        end                        
                    end
                end
            end
                
            if P.cellForMultipleTrials
                GDF = {};
                for trial = 1:length(P.trialIDs)
                    for unit = 1:length(P.unitIDs)
                        if ~isempty(spiketrain_container{trial,unit})
                            GDF{trial, unit} = spiketrain_container{trial,unit}(:,2);
                        else
                            GDF{trial, unit} = [];
                        end
                    end
                end
            else
                GDF = [];
                for trial = 1:length(P.trialIDs)
                    tgdf = [];
                    for unit = 1:length(P.unitIDs)
                        tgdf = [tgdf; spiketrain_container{trial, unit}];
                    end
                    if ~P.alignOnEvent && ~isempty(tgdf)
                        tgdf(:,2) = tgdf(:,2) + (trial-1)*15*32000;
                    end
                    if ~isempty(tgdf)
                        tgdf = sortrows(tgdf,2);
                        GDF = [GDF ; tgdf];
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        function R = getGDF_(self, P)
            gdf = true; cut = false;
            
            if ~isempty(P.trialIDs) && ~isempty(P.unitIDs)
                constraint_str = [];
                if ~isempty(P.from); constraint_str = [constraint_str ' AND s.sample>=' num2str(P.from)]; end
                if ~isempty(P.to);   constraint_str = [constraint_str ' AND s.sample<=' num2str(P.to)]; end 
                qstr = ['SELECT u.index, s.sample FROM spike as s JOIN segment as t ON (s.segment_id = t.segment_id) JOIN unit as u ON (u.unit_id = s.unit_id) WHERE s.segment_id IN (' util.dlmstring(P.trialIDs) ') AND s.unit_id IN (' util.dlmstring(P.unitIDs) ') ' constraint_str ' ORDER BY t.index, s.sample, u.index'];
            elseif (~isempty(P.trialIDs) && ~isempty(P.analysisIDs))
                [constraint_str] = self.prepareGDFStatements(gdf, cut, P);
                select_str = 'SELECT u.index, s.sample FROM analysis as a JOIN unit as u ON (u.analysis = a.analysis_id) JOIN spike as s ON (s.unit_id = u.unit_id)';
                qstr = [select_str ' WHERE s.segment_id IN (' util.dlmstring(P.trialIDs) ') AND a.analysis_id IN (' util.dlmstring(P.analysisIDs) ') ' constraint_str ' ORDER BY s.sample ASC'];
            else
                error('not implemented');             
            end
            R = cell2mat(self.query(qstr, {'Double', 'Double'}));          
        end      
        
        %------------------------------------------------------------------
        function EGDF = getEGDF(self, varargin)
            P = self.perpareInputs([], varargin{:});
            % If multiple trials are selected the sql statement down there
            % gets insanely slow. fix this or live with this workaround:
            import matlabfilecentral.progressbar.progressbar
            progressbar(0, 0)         % Initialize/reset two bars
            progressbar('Trials', 'Units');
            EGDF = [];
            P_ = P;            
            for trial = 1:length(P.trialIDs)
                for unit = 1:length(P.unitIDs)
                    P_.segment_idIDs = P.trialIDs(trial);
                    P_.unitIDs  = P.unitIDs(unit);
                    
                    pHash = self.P2Hash(P_, 'unitIDs', 'trialIDs', 'analysisIDs', 'cutleft', 'Tf', 'from', 'to', 'channelIDs', 'channelNRs', 'unitIDs', 'tetrodeIDs');
                    egdf = self.EGDFBuffer.call(pHash, P_);    
                    EGDF = [EGDF; egdf];
                    progressbar([], unit/length(P.unitIDs));
                end
                progressbar(trial/length(P.trialIDs), []);
            end
        end
        %------------------------------------------------------------------
        function R = getEGDF_(self, P)
            [constraint_str] = self.prepareGDFStatements(true, true, P);
            
            if (~isempty(P.trialIDs) && ~isempty(P.analysisIDs))
                select_str = 'SELECT u.index, s.sample, c.index, td.data[s.sample-%d:s.sample+%d] FROM analysis as a JOIN tetrode as te ON (a.channel_group_id = te.experiment_id) JOIN channel as c ON (c.tetrode = te.experiment_id) JOIN trialdata as td ON (td.channel = c.channel_id) JOIN unit as u ON (u.analysis = a.analysis_id) JOIN spike as s ON (s.unit_id = u.unit_id AND s.segment_id = td.segment_id)';
                qstr = [select_str ' WHERE s.sample-%d>0 AND td.segment_id IN (' util.dlmstring(P.trialIDs) ') ' constraint_str ' ORDER BY u.index, s.sample, c.index ASC'];
                qstr = sprintf(qstr, P.cutleft, P.Tf-P.cutleft-1, P.cutleft);                
            else
                error('not implemented');             
            end
            R = self.query(qstr);
            % deal with partially cut spikes at the beginning and end of
            % the timeseries
            if isempty(R)
                R = [];
                return
            end
            nC = length(unique(cell2mat(R(:,3)))); 

            
            for i=size(R,1):-1:1
                if length(R{i,4}) < length(R{1,4})
                    R(i,:) = [];
                end
            end
            R = [cell2mat(R(:,1:3)) cell2mat(R(:,4)')'];  
            if size(R,2) == 1
                R = R';
            end
            rr = R;
            R = [];
            count = 1;
            for i=1:nC:size(rr,1)
                R(count,:) = [rr(i,1:2) mysort.util.m2v(rr(i:i+nC-1, 4:end))];
                count = count +1;
            end
        end
        %------------------------------------------------------------------
        function R = getCutSpikes(self, varargin)
            R = self.getEGDF(varargin{:});
            R = R(:,3:end);
        end
        
        %------------------------------------------------------------------
        function [N pos] = getNoiseSamples(self, varargin)
            P.Len = 100;
            P.Number = 10;
            P.minSpikeDistance = 50;
            P.gdf = [];
            P = self.perpareInputs(P, varargin{:});
            
            nC = length(P.channelIDs);
            if isempty(P.gdf)
                [P.gdf null_trials] = self.getGDF(varargin{:});
            end
               
            if isempty(P.to)
                R = self.getTrialLength(varargin{:});
                trial_len = R(2);
            else
                trial_len = P.to;
            end
            spike_times = [P.from; sort(P.gdf(:,2)); trial_len];
            
            % compute interspike intervals to find periods of noise
            isi = [diff(spike_times) spike_times(1:end-1)];
                        
            % remove isi that are too short
            isi(isi(:,1) <= 2*P.minSpikeDistance+P.Len, :) = [];
            assert(~isempty(isi), 'There are no noise periods of that length!');
            
            % Check if enough noise is left
            cumsum_isi = cumsum(isi(:,1));
            noiseTotalLen = cumsum_isi(end);
            assert(P.Number*P.Len < 3*noiseTotalLen, 'The total amount of noise is covered three times by requested samples!');
            
            N = zeros(P.Number, nC*P.Len);
            pos = zeros(P.Number, 1);

            if P.Number > 300
                XX = self.getData(varargin{:}, 'from', 1, 'to', trial_len);
            end
            for k=1:P.Number
                % choose randomly a noise epoch with probability given by
                % noise epoch length
                r = rand()*(noiseTotalLen-1);
                epoch_idx = find(cumsum_isi>r, 1);
                
                % choose randomly a position inside that noise epoch
                % which has enough distance to the next spikes
                start_epoch = isi(epoch_idx,2) +   P.minSpikeDistance;
                len_epoch   = isi(epoch_idx,1) - 2*P.minSpikeDistance;
                assert(len_epoch > 0, 'Epoch Len is zero or negative!');
                offset_inside_epoch = round(rand()*(len_epoch-P.Len));
                assert(offset_inside_epoch >= 0, 'offset_inside_epoch is negative!');
                pos(k) = start_epoch + offset_inside_epoch;
                
                violated_spikes = find( (P.gdf(:,2) >= pos(k)) & ...
                                        (P.gdf(:,2) <= pos(k)+P.Len-1) );
                assert(isempty(violated_spikes), 'Spikes are inside epoch!');
                % cut the epoch
                if P.Number > 300
                    X = XX(:, pos(k):pos(k)+P.Len-1);
                else
                    X = self.getData(varargin{:}, 'from', pos(k), 'to', pos(k)+P.Len-1);
                end
                                
                N(k,:) = mysort.util.m2v(X);
            end
        end
        
        %------------------------------------------------------------------
        function [constraint_str] = prepareGDFStatements(self, gdf, cut, P)
            if cut
                assert(~isempty(P.Tf) && ~isempty(P.cutleft), 'cutleft and Tf must be specified!');
            end
            constraint_str = [];
            if ~isempty(P.from); constraint_str = [constraint_str ' AND s.sample>=' num2str(P.from)]; end
            if ~isempty(P.to);   constraint_str = [constraint_str ' AND s.sample<=' num2str(P.to)]; end 

            if cut && ~isempty(P.channelIDs)
                constraint_str = [constraint_str ' AND c.channel_id IN (' util.dlmstring(P.channelIDs) ')'];
            elseif cut && ~isempty(P.channelNRs)
                constraint_str = [constraint_str ' AND c.index IN (' util.dlmstring(P.channelNRs) ')'];
            end
            if ~isempty(P.unitIDs)
                constraint_str = [constraint_str ' AND u.unit_id IN (' util.dlmstring(P.unitIDs) ')'];
            end
%             if ~isempty(P.tetrodeIDs)
%                 constraint_str = [constraint_str ' AND te.experiment_id IN (' util.dlmstring(P.tetrodeIDs) ')'];
%             end
        end
        %------------------------------------------------------------------
        function [T, algoIDs] = getTemplatesFromCutSpikes(self, varargin)
            P.method = 'mean';
            P = self.perpareInputs(P, varargin{:});
            
            rr = self.getEGDF('otherP', P);
            algoIDs = unique(rr(:,1));
            T = zeros(length(algoIDs), size(rr,2)-2);
            if ~isempty(P.method) && strcmp(P.method, 'median')
                fun = @median;
            else
                fun = @mean;
            end
            for i=1:length(algoIDs)
                T(i,:) = fun(rr(rr(:,1)==algoIDs(i), 3:end),1);
            end
        end
        %------------------------------------------------------------------
        function [T trialIDs trialidx unitIDs algoIDs cutleft snr] = getTemplatesConcat(self, varargin)
            P = self.perpareInputs([], varargin{:});
            if ~isempty(P.channelIDs); channel_constraint = ['AND c.channel_id IN (' util.dlmstring(P.channelIDs) ') ']; else channel_constraint=[]; end
            qstr = ['SELECT wfs.tid, wfs.tidx, wfs.uid, wfs.algoid, wfs.cutleft, wfs.snr, wfs.index, wfs.data FROM (SELECT t.segment_id as tid, t.index as tidx, u.unit_id as uid, u.index as algoid, c.index, wf.data as data, tmp.cutleft as cutleft, tmp.snr as snr FROM unit AS u JOIN template AS tmp ON tmp.unit_id = u.unit_id JOIN segment AS t ON t.segment_id = tmp.segment_id JOIN waveform as wf ON wf.template_id = tmp.template_id join channel as c on c.channel_id = wf.channel_id where u.unit_id IN (' util.dlmstring(P.unitIDs) ') and t.segment_id IN (' util.dlmstring(P.trialIDs) ') ' channel_constraint 'ORDER BY t.index, u.index, u.unit_id, c.index asc) as wfs '];
            R = self.query(qstr, {'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray'});
            if isempty(R);
                T =[]; trialIDs =[]; trialidx =[]; unitIDs =[]; algoIDs =[]; cutleft =[]; snr =[];
                return
            end
            nC = length(unique(cell2mat(R(:,7))));
            T = cell2mat(R(:,8)')';
            if nC == 0
                nC = size(T,2)/P.Tf;
            end
            R = R(1:nC:end,:);
            trialIDs   = cell2mat(R(:,1));
            trialidx   = cell2mat(R(:,2));
            unitIDs    = cell2mat(R(:,3));
            algoIDs    = cell2mat(R(:,4));
            cutleft    = cell2mat(R(:,5));
            snr        = cell2mat(R(:,6));            

            if size(T,2) > P.Tf
                T = T(:,1:P.Tf);
            elseif size(T,2) < P.Tf
                T = [T zeros(size(T,1), P.Tf-size(T,2))];
            end    
            T = mysort.util.m2v(T, nC);
           
        end
        %------------------------------------------------------------------
        function [T trialIDs trialidx unitIDs algoIDs channelIDs channelNRs cutleft snr] = getTemplatesChannel(self, varargin)
            P.minmax = 0;
            P = self.perpareInputs(P, varargin{:});
            channelconstr = [];
            if ~isempty(P.channelIDs)
                channelconstr = [' AND c.channel_id IN (' util.dlmstring(P.channelIDs) ') '];
            elseif ~isempty(P.channelNRs)
                channelconstr = [' AND c.index IN (' util.dlmstring(P.channelNRs) ') '];
            end
            
            if P.minmax
                wavdatastr = 'array_minmax(wav.data)';
            else
                wavdatastr = 'wav.data';
            end
            
            if ~isempty(P.unitIDs) && ~isempty(P.trialIDs)
                qstr = ['SELECT t.segment_id, t.index, u.unit_id, u.index, c.channel_id, c.index, tmp.cutleft, tmp.snr, ' wavdatastr ' FROM unit as u JOIN template as tmp ON (tmp.unit = u.unit_id) JOIN waveform as wav ON (tmp.template_id = wav.template) JOIN channel as c ON (c.channel_id = wav.channel) JOIN segment as t ON (tmp.segment_id = t.segment_id) WHERE t.segment_id in (' util.dlmstring(P.trialIDs) ')  AND u.unit_id IN (' util.dlmstring(P.unitIDs) ') ' channelconstr ' ORDER BY ( t.index, u.unit_id, c.index)'];
            else
                error('not implemented');
            end         
            R = self.query(qstr, {'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray'});
            trialIDs   = cell2mat(R(:,1));
            trialidx   = cell2mat(R(:,2));
            unitIDs    = cell2mat(R(:,3));
            algoIDs    = cell2mat(R(:,4));
            channelIDs = cell2mat(R(:,5));
            channelNRs = cell2mat(R(:,6));
            cutleft    = cell2mat(R(:,7));
            snr        = cell2mat(R(:,8));
            T = cell2mat(R(:,9)')';
        end
        %------------------------------------------------------------------
        function R = getAnalysis(self, varargin)
            P = self.perpareInputs([], varargin{:});
            constraints = {};
            if ~isempty(P.analysisIDs)
                constraints = [constraints ['a.analysis_id IN (' util.dlmstring(P.analysisIDs) ')']];
            end
%             if ~isempty(P.blockIDs)
%                 constraints = [constraints ['a.block IN (' util.dlmstring(P.blockIDs) ')']];
%             end
            if ~isempty(P.experimentIDs)
                constraints = [constraints ['cg.experiment_id IN (' util.dlmstring(P.experimentIDs) ')']];
            end
            if ~isempty(P.segment_idIDX)
                constraints = [constraints ['a.idx_start <= ' num2str(max(P.segment_idIDX)) ' AND a.idx_end >= ' num2str(min(P.segment_idIDX)) ' ']];
            end
            if ~isempty(P.tetrodeIDs)
                constraints = [constraints ['a.channel_group_id IN (' util.dlmstring(P.tetrodeIDs) ')']];
            end            
            qstr = 'SELECT a.analysis_id, aa.value, a.created_on, a.kind_code FROM analysis AS a JOIN analysis_attribute AS aa ON (a.analysis_id=aa.analysis_id) JOIN channel_group AS cg ON (cg.channel_group_id = a.channel_group_id)';
            if ~isempty(constraints)
                qstr = [qstr 'WHERE aa.name IN (''name'') AND ' util.dlmstring(constraints, ' AND ', '%s') ' ORDER BY a.kind_code, a.created_on, aa.value'];
            end
            R = self.query(qstr, {'Double', 'Str', 'Str', 'Str'});            
        end
        %------------------------------------------------------------------
        function R = getAnalysisDetails(self, analysisIDs)
            R = {};
            if isempty(analysisIDs); return; end
            qstr = ['SELECT * FROM analysis WHERE id IN (' util.dlmstring(analysisIDs) ')'];            
            R = self.query(qstr, {'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Str', 'Str', 'Str', 'Str', 'Str', 'Str'});        
        end
        
        %------------------------------------------------------------------
        function R = getUnits(self, varargin)
            P = self.perpareInputs([], varargin{:});
            if ~isempty(P.analysisIDs)
                qstr = ['SELECT u.unit_id, u.index, u.kind_code FROM unit AS u WHERE u.analysis_id IN (' util.dlmstring(P.analysisIDs) ') ORDER BY u.index, u.unit_id ASC'];
            else
                error('not implemented');
            end
            R = self.query(qstr, {'Double', 'Double', 'Str'});            
        end
        
        %------------------------------------------------------------------
        function R = getEventTypes(self)
            qstr = 'SELECT kind_code FROM event GROUP BY kind_code';
            R = self.query(qstr, {'Double','Str'});  
        end
        %------------------------------------------------------------------
        function R = getEvents(self, varargin)
            %P.returnEmptyTrials = 1;
            P = self.perpareInputs([], varargin{:});
            
            subselectstr = {};
            if ~isempty(P.to)
                subselectstr = [subselectstr ['es.sample<=' num2str(P.to)]];
            end
            if ~isempty(P.from)
                subselectstr = [subselectstr ['es.sample>=' num2str(P.from)]];
            end
            if ~isempty(P.eventIDs)
                subselectstr = [subselectstr ['e.experiment_id IN (' util.dlmstring(P.eventIDs) ')']];
            end
            subselectstr = util.dlmstring(subselectstr, ' AND ', '%s');
            
            if length(P.trialIDs) == 1
                if ~isempty(subselectstr); subselectstr = [subselectstr ' AND ']; end
                qstr = ['SELECT es.id, es.sample, e.experiment_id, e.name FROM event as e JOIN eventsample AS es ON (es.event = e.experiment_id) WHERE ' subselectstr ' es.segment_id =' num2str(P.trialIDs)];
                R = self.query(qstr, {'Double', 'Double', 'Double', 'Str'});   
            elseif ~isempty(P.trialIDs) && length(P.eventIDs) == 1 
                qstr = ['SELECT es.id, es.sample, e.experiment_id, e.name, t.segment_id FROM event as e JOIN eventsample AS es ON (es.event = e.experiment_id) JOIN trial AS t ON (t.segment_id = es.segment_id) WHERE t.segment_id IN (' util.dlmstring(P.trialIDs) ') AND ' subselectstr ' ORDER BY t.index'];
                R = self.query(qstr, {'Double', 'Double', 'Double', 'Str', 'Double'});   
%             elseif ~isempty(P.trialIDs) && length(P.eventIDs) == 1 && P.returnEmptyTrials
%                 qstr = ['SELECT es.id, es.sample, e.experiment_id, e.name FROM segment as t LEFT OUTER JOIN eventsample AS es ON (es.segment_id = t.segment_id) LEFT OUTER JOIN event AS e ON (es.event = e.experiment_id) WHERE t.segment_id IN (' util.dlmstring(P.trialIDs) ') AND (' subselectstr ' OR es.segment_id IS NULL) ORDER BY t.index'];
%             elseif ~isempty(P.trialIDs) && length(P.eventIDs) == 1 && ~P.returnEmptyTrials
%                 qstr = ['SELECT es.id, es.sample, e.experiment_id, e.name FROM event as e JOIN eventsample AS es ON (es.event = e.experiment_id) JOIN trial AS t ON (t.segment_id = es.segment_id) WHERE t.segment_id IN (' util.dlmstring(P.trialIDs) ') AND ' subselectstr ' ORDER BY t.index'];
            else
                error('not implemented');
            end
                     
        end
        %------------------------------------------------------------------
        function E = getEventsStruct(self, varargin)
            events = self.getEvents(varargin{:});
            E = [];
            for i=1:size(events,1)
                E(i).name    = events{i,4};
                E(i).sample  = events{i,2};
                E(i).channel = events{i,3};
            end
        end
        %------------------------------------------------------------------
        function [X T]= getData(self, varargin)
            P.downsampleTo = [];
            P = self.perpareInputs(P, varargin{:});
            T = []; subselectstr=[];
            if ~isempty(P.to) && ~isempty(P.from)
                subselectstr = ['[' num2str(P.from) ':' num2str(P.to) ']'];
            elseif ~isempty(P.to)
                subselectstr = ['[1:' num2str(P.to) ']'];
            elseif ~isempty(P.from)
%                 subselectstr = ['[' num2str(P.from) ':-1]'];
                error('At the moment it is not possible to only provide "from" without "to" !');
            end
            
            if ~isempty(P.downsampleTo)
                subselectstr = ['downsample_data_for_plotting(td.data' subselectstr ', ' num2str(P.downsampleTo) ')']; 
            else
                subselectstr = ['td.data' subselectstr];
            end
            X = [];
            for trialIdx=1:length(P.trialIDs)
                if ~isempty(P.trialIDs) && ~isempty(P.tetrodeIDs)
                    qstr = ['SELECT c.index, ' subselectstr ' FROM analog_signal as td JOIN channel as c ON (c.channel_id = td.channel_id) JOIN channel_group as te ON (c.channel_group_id=te.channel_group_id) WHERE td.segment_id = %d AND te.channel_group_id = %d ORDER BY c.index ASC'];
                    qstr = sprintf(qstr, P.trialIDs(trialIdx), P.tetrodeIDs); 
                else
                    error('not implemented');
                end
                R = self.query(qstr, {'Double', 'DoubleArray'});
                XX = [];
                for i=1:size(R,1)
                    if isempty(R{i,2})
                        warning('empty result set, maybe from-to boundaries not correct?');
                        return
                    end
                    XX(i,:) = R{i,2};
                end

                if ~isempty(P.downsampleTo)
                    T = XX(:,end/2+1:end);
                    XX = XX(:,1:end/2);
                end  
                X = [X XX];
            end
        end
                
        %------------------------------------------------------------------
%         function X = getAmplitudeDrift(self, varargin)
%             P = self.perpareInputs([], varargin{:});
%             if ~isempty(P.unitIDs) && ~isempty(P.trialIDs)
%                 for i=1:size(X,1)
%                     
%                 end
%             else
%                 error('not implemented');
%             end
%         end   
        %------------------------------------------------------------------
        function anaID = insertSorting(self, monkeyName, expName, blockName, tetrodeNumber, trialIdxList,...
                algoName, algoNeuronIDs, gdfList, templateList, cutleft, covMatrixList)
            % Inserts a complete sorting into the DB
            %
            % Inputs:
            %   monkeyName     - the name of the monkey, "Julia", "Louis"
            %   expName        - the name of the experiment, "L011"
            %   blockName      - name of the block to which the trials
            %                    belong
            %   tetrodeNumber  - number of the tetrode (1-16)
            %   trialIdxList   - Array containing the trial indices, 1:20
            %   algoName       - the name of the algorithm that was used to
            %                    create this sorting
            %   algoNeuronIDs  - array with the algo ids of the templates
            %
            %                   WARNING!!!
            %   The algoNeuron ID in the DB starts at 0 !! 
            %  
            %
            %   gdfList        - cell array containing a gdf per trial in
            %                    trialIdxList. The IDs of neurons in the
            %                    gdf are the algo_ids in the DB!
            %   templateList   - cell array containing the templates for
            %                    each trial in trialIdxList
            %                    OR
            %                    tensor containing the templates for the
            %                    whole analysis, if they dont change over
            %                    the trials. dimensions are
            %                    samples x channels x neurons
            %   cutleft
            %   covMatrixList -  cell array containing the covariance
            %                    functions for all pairs of channels for
            %                    all trials in trialIdxList
            %                    OR
            %                    a matrix containing the xcorr function in
            %                    the way matlab stores them (see xcorr
            %                    function)
            %
            % Output:
            %   anaId  - the analysis ID of the sorting that was just
            %            inserted.
            assert(length(trialIdxList) == length(gdfList), 'There must be one trialIdx per gdf!');
            if iscell(templateList)
                T = templateList{1};
            else
                T = templateList;
            end
            [Tf nC nT] = size(T);
            assert(length(algoNeuronIDs) == nT, 'there must be one template per neuron!');
            
            R = self.getExperiments('monkeys', monkeyName, 'experiments', expName);
            assert(~isempty(R), 'Could not find experiment!');
            expId = R{1};
            R = self.getTrials('experimentIDs', expId);
            assert(~isempty(R), 'No trials for this experiment!');
            R = cell2mat(R);
            trialIDs = R(trialIdxList,1);
            R = self.getTetrodes('experimentIDs', expId, 'tetrodeNRs', tetrodeNumber);
            assert(~isempty(R), 'Tetrode number not found!');
            tetrodeID = R{1,1};
            R = self.getBlocks('experimentIDs', expId, 'blocks', blockName);
            assert(~isempty(R), 'Block not found!');
            blockID = R{1,1};
            % Get new analysis ID
            anaID = self.insertNewAnalysis(expId, blockID, tetrodeID, min(trialIdxList), max(trialIdxList), algoName);

            % Insert new neurons and the spike trains
            unitIDsAlgoIDs = insertSortingAndUnits(self, anaID, trialIDs, gdfList, algoNeuronIDs);
            
            R = self.getChannels('tetrodeIDs', tetrodeID);
            channelIDs = cell2mat(R(:,1)); 
            
            % add templates to units
            self.insertTemplates(unitIDsAlgoIDs, trialIDs, templateList, cutleft, channelIDs);
                
            % add covariance matrices            
            self.insertCovariance(anaID, trialIDs, channelIDs, 'noise', covMatrixList);
        end
        
        %------------------------------------------------------------------
        function anaId = insertNewAnalysis(self, expId, blockId, tetrodeId, trialstart, trialend, algo_main)
            str = 'INSERT INTO analysis (expId, block, tetrode, trialidxstart, trialidxend, algorithm, kind, status) VALUES ';
            str = [str sprintf('(%d, %d, %d,%d,%d,''%s'', ''SORT'', ''auto_filling'')', ...
                expId, blockId, tetrodeId, trialstart, trialend, algo_main)];
%             disp(str)
            self.DBH.execute(str);
            R = self.query('SELECT currval(pg_get_serial_sequence(''analysis'', ''id''))');
            anaId = R{1,1};
        end
        %------------------------------------------------------------------
        function unitIDsAlgoIDs = insertSortingAndUnits(self, anaId, trialIds, gdfList, algoNeuronIDs)
            assert(length(trialIds) == length(gdfList), 'There must be one trialID per gdf!');
            
            existingUnitIDsAlgoIDs = self.getUnits('analysisIDs', anaId);
            if ~isempty(existingUnitIDsAlgoIDs)
                existingUnitIDsAlgoIDs = cell2mat(existingUnitIDsAlgoIDs(:,1:2));
            else
                existingUnitIDsAlgoIDs = [];
            end
            % First, make sure that all neurons exist in the DB
            if isempty(existingUnitIDsAlgoIDs)
                missingIDs = algoNeuronIDs;
            else
                missingIDs = setdiff(algoNeuronIDs, existingUnitIDsAlgoIDs(:,2));
            end
            if ~isempty(missingIDs)
                existingUnitIDsAlgoIDs = [existingUnitIDsAlgoIDs; self.insertNewUnits(anaId, missingIDs, 'single')];
            end
            % sort the unitID vs algoID mapping to that the algo IDs have
            % the same order as in algoNeuronIDs
            unitIDsAlgoIDs = zeros(size(existingUnitIDsAlgoIDs));
            for i=1:length(algoNeuronIDs)
                euidx = find(existingUnitIDsAlgoIDs(:,2) == algoNeuronIDs(i));
                unitIDsAlgoIDs(i,:) = existingUnitIDsAlgoIDs(euidx,:);
            end
            for i=1:length(trialIds)
                % replace the algoIds in the gdf with the database unit ids
                gdf = gdfList{i};
                ugdf = gdf; % make a copy to not overwrite the ids all the time
                for u = 1:size(existingUnitIDsAlgoIDs)
                    ugdf(gdf(:,1)==existingUnitIDsAlgoIDs(u,2),1) = existingUnitIDsAlgoIDs(u,1);
                end
                
                % Now insert spike train
                self.insertSpikeTrain(trialIds(i), ugdf);
            end
        end                
        %------------------------------------------------------------------
        function unitIds = insertNewUnits(self, anaId, algoIds, kind)
            str = 'INSERT INTO unit (analysis, algoid, kind) VALUES (%d, %d, ''%s'')';
            unitIds = zeros(length(algoIds),2);
            for i = 1:length(algoIds)
                q_str = sprintf(str, anaId, algoIds(i), kind);
                disp(q_str)
                self.DBH.execute(q_str);
                R = self.query('SELECT currval(''unit_id_seq'')');
                unitIds(i,:) = [R{1,1} algoIds(i)];
            end
        end
        %------------------------------------------------------------------
        function insertSpikeTrain(self, trialId, gdf)
            str = 'INSERT INTO spike (trial, unit, sample) VALUES ';
            for i=1:size(gdf,1)
                str = [str sprintf('(%d, %d, %d),', trialId, gdf(i,1), gdf(i,2))];
            end
            str = str(1:end-1);
            self.DBH.execute(str);
        end                
        %------------------------------------------------------------------
        function tIds = insertTemplates(self, unitIDsAlgoIDs, trialIDs, T, cutleft, channelIDs)
            if iscell(T)
                assert(length(T) == length(trialIDs), 'There must be as many templates as trials!');
                nT = size(T{1},3);
                nC = size(T{1},2);
                Tf = size(T{1},1);
            else
                nT = size(T,3);
                nC = size(T,2);
                Tf = size(T,1);
            end
            assert(nC == length(channelIDs), 'nC does not match with channels!');
            tIDs = zeros(length(trialIDs), nT);
            template_str = 'INSERT INTO template (unit, trial, cutleft) VALUES (%d, %d, %d)';
            wf_str = 'INSERT INTO waveform (template, channel, data) VALUES ';
            for c=1:nC
                wf_str = [wf_str sprintf('(%%d, %d, ''{%%%%s}''),', channelIDs(c))];
            end
            wf_str = wf_str(1:end-1);
            for u=1:size(unitIDsAlgoIDs,1)
                myUID = unitIDsAlgoIDs(u,1);
                myAID = unitIDsAlgoIDs(u,2);
                if ~iscell(T)
                    % The unit has only one template for all trials, insert
                    % always the same one
                    data_str = {};
                    for c = 1:nC
                        data_str{c} = sprintf('%.6f,', squeeze(T(:,c,u)));
                        data_str{c} = data_str{c}(1:end-1);
                    end    
                end
                for i=1:length(trialIDs)
                    temp_q_str = sprintf(template_str, myUID, trialIDs(i), cutleft);
                    self.DBH.execute(temp_q_str);
                    R = self.query('SELECT currval(''template_id_seq'')');
                    tempID = R{1,1};
                    tempIDChannelCell = mat2cell(repmat(tempID,1,nC),1,ones(1,nC));
                    wf_t_str = sprintf(wf_str, tempIDChannelCell{:});
                    if iscell(T) 
                        % The unit has different templates for every trial,
                        % build everytime a new data str
                        data_str = {};
                        for c = 1:nC
                            data_str{c} = sprintf('%.6f,', squeeze(T{i}(:,c,u)));
                        end
                    end
                    wf_q_str = sprintf(wf_t_str, data_str{:});
                    self.execute(wf_q_str);
                    tIDs(i,u) = tempID;
                end
            end
        end                
        %------------------------------------------------------------------
        function covId = insertCovariance(self, anaId, trialIDs, channelIDs, kind, CList)
            if ~iscell(CList)
                CList = {CList};
            end
            for i=1:length(CList)
                t_id = trialIDs(i);
                C = CList{i};
                cov_istr = 'INSERT INTO covariance (analysis, trial, kind) VALUES (%d, %d, ''%s'')';
                q_cov_istr = sprintf(cov_istr, anaId, t_id, kind);
                self.DBH.execute(q_cov_istr);

                req_str = 'SELECT currval(''covariance_id_seq'');';
                R = self.query(req_str);
                covId = R{1,1};

                xcorr_istr = 'INSERT INTO xcorr (covariance, channel1, channel2, data) VALUES (%d, %d, %d, ''{%s}'')';
                nC = length(channelIDs);
                count = 1;
                for c1 = 1:nC
                    for c2 = 1:nC
                        if c2>=c1
                            dstr = sprintf('%.6f,', C(:,count));
                            dstr = dstr(1:end-1);
                            q_xcorr_istr = sprintf(xcorr_istr, covId, channelIDs(c1), channelIDs(c2), dstr);
                            self.DBH.execute(q_xcorr_istr);
                        end
                        count = count+1;
                    end
                end
            end
        end                

        %------------------------------------------------------------------
        function P = perpareInputs(self, P, varargin)
            P.monkeys = [];
            P.experiments = [];
            P.experimentIDs = [];
            P.blocks = [];
            P.blockIDs = [];
            P.segment_idIDX = [];
            P.trialIDs = [];
            P.tetrodeIDs = [];
            P.tetrodeNRs = [];
            P.channelIDs = [];
            P.channelNRs = [];
            P.analysisIDs = [];
            P.analysisKind = [];
            P.unitIDs = [];
            P.from = [];
            P.to   = [];
            P.waveformIDs = [];
            P.templateIDs = [];
            P.cutleft = [];
            P.Tf = [];
            P.events = [];
            P.eventIDs = [];
            
            P = mysort.util.parseInputs(P, varargin, 'merge');
        end
        %------------------------------------------------------------------
        function hash = P2Hash(self, P, varargin)
            hash = [];
            for i=1:length(varargin)
                p = P.(varargin{i});
                if size(p,1)>0
                    p = p';
                end
                hash = [hash '_' num2str(p)];
            end
            hash = int64(java.lang.String(hash).hashCode);
        end
    end
end