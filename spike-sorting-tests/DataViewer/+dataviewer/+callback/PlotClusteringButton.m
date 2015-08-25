function handles = PlotClusteringButton(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    
    if length(P.analysisIDs) > 1
        handles.warning_func('Only one analysis might be selected at a time!');
        return
    end
    C_cell = handles.DH.getCovarianceMatrixCellStruct('analysisIDs', P.analysisIDs, 'trialIDs', P.trialIDs(1));
    if isempty(C_cell)
        handles.warning_func('No Noise Covariance at first selected trial!');
        return
    end
    if C_cell{1,1}.Tf < P.Tf
        handles.warning_func(['You selected a Tf (' num2str(P.Tf) ') that is bigger than that of the stored noise covariance matrix (' num2str(C_cell{1,1}.Tf) ')!']);
        return
    end

    N = mysort.util.NoiseEstimator(C_cell);
    C = N.getNoiseCovarianceMatrix(P.Tf, PP.channel.idx);
    C = C + C(1,1)*eye(size(C));
%     [T, trialidx, trialIDs, unitIDs, algoID, cutleft] = handles.DH.getTemplatesConcat('otherP', P);
%     if isempty(T) || size(T,1) < length(P.unitIDs)
%         handles.warning_func('At least one selected Units has no Template!!!');
%         return
%     end
%     T(algoID+1,:) = T;

    
    X = handles.DH.getEGDF('otherP', P);
    IDs = unique(X(:,1));
    for t=1:length(IDs)
        T(t,:) = mean(X(X(:,1)==IDs(t),3:end),1);
    end
    
    try
        U = chol(C);
    catch
        handles.warning_func('Could not compute Choleskyfactorization of C!');
        return
    end
    iU = inv(U);
    pwSpikesX   = X(:,3:end)*iU; 
    pwTemplates = T*iU;
    
    fetX = mysort.util.dimReductionPCA([pwSpikesX; pwTemplates], 6);
    fetSpikes    = fetX(1:end-size(T,1),:);
    fetTemplates = fetX(end-size(T,1)+1:end,:);
    mysort.plot.clustering(fetSpikes, X(:,1), fetTemplates, IDs);    
    mysort.plot.figureName('Clustering'); 
    mysort.plot.figureTitle(P.figureTitle); 
    handles.new_figure_handles = gcf;
    