function handles = acceptTemplate(handles)
    G = handles.GUI;
    G.acceptCurrentTemplate();
    G.update(handles);
    
        
        


%         handles.DATA = hdmeagui.data.calcSpikeAmplitudes(handles.DATA);
%         handles.VIEW = hdmeagui.view.plotChannelAmp(handles);
%         handles.DATA = hdmeagui.data.update(handles);
%         handles.VIEW = hdmeagui.view.update(handles);   
%     end