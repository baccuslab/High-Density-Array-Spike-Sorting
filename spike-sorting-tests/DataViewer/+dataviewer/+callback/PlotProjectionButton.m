function handles = PlotProjectionButton(hObject, eventdata, handles)
    [P PP] = dataviewer.util.getHandleValues(handles);
    nC = 4;
    C_cell = handles.DH.getCovarianceMatrixCellStruct('analysisIDs', P.analysisIDs, 'trialIDs', P.trialIDs(1));
    N = mysort.util.NoiseEstimator(C_cell);
    C = N.getNoiseCovarianceMatrix(P.Tf, PP.channel.idx);
    C = C + C(1,1)*eye(size(C));
% 
%     [T, trialidx, trialIDs, unitIDs, algoID, cutleft] = handles.DH.getTemplatesConcat('otherP', P);
%     if isempty(T) || size(T,1) < length(P.unitIDs)
%         handles.warning_func('At least one selected Units has no Template!!!');
%         return
%     end
%     T(algoID+1,:) = T;
   
    
    X = handles.DH.getEGDF('otherP', P);
    
    mysort.plot.clusterProjection(X(:,3:end), X(:,1), [], C);
    mysort.plot.figureName('Projection Plot'); 
    mysort.plot.figureTitle(P.figureTitle); 
    handles.new_figure_handles = gcf;
 