function hidens_write_neuroplacement(configPath, fname, npos, n_size, varargin)

    xsp = n_size(1);
    ysp = n_size(2);
    
    P.folder = 'matlab_specs';

    % larger size for wider configs
    P.size=round(max(xsp,ysp)*20);
    
    P.elcnt=nan;
    P.multiloc=0;
    P = mysort.util.parseInputs(P, varargin);
    fullspecsfolder = fullfile(configPath, P.folder);
    if not(exist(configPath, 'dir'))
        error('Configs directory does not exist');
    end
    if not(exist(fullspecsfolder, 'dir'))
        mkdir(fullspecsfolder);
    end
    full_fname = fullfile(fullspecsfolder, fname);


    %export to file
    fid=fopen(full_fname, 'w');
    for i=1:length(npos)
        dx = P.size;  %size
        if isfield(npos{i}, 'dx')
            dx=npos{i}.dx;
        end
        
        dy = dx;
        if isfield(npos{i}, 'dy')
            dy = npos{i}.dy;
        end
        
        if P.multiloc
            fprintf(fid, '*');
        end
        
        fprintf(fid, 'Neuron %s: %d/%d, %d/%d', npos{i}.label,...
            round(npos{i}.x), round(npos{i}.y), round(dx), round(dy));
        
        if isfield(npos{i}, 'c_sr')
            fprintf(fid, ', sr%d', npos{i}.c_sr);
        end
        if isfield(npos{i}, 'cost')
            fprintf(fid, ', c%f', npos{i}.cost);
        end
        
        elcnt = P.elcnt;
        if isfield(npos{i}, 'elcnt')
            elcnt=npos{i}.elcnt;
        end
        if ~isnan(elcnt)
            fprintf(fid, ', el%d', elcnt);
        end
        if npos{i}.isStim
            fprintf(fid, ', stim');
        end        
        fprintf(fid, '\n');
    end
    fclose(fid);












