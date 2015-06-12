clear all
close all
clear all
fname = 'C:\LocalData\Michele\Trace_id858_2012-03-06T09_50_14_0to9.h5';
gui = struct();
gui.Window = mysort.plot.figure('w', 800, 'h', 500, 'name', 'H5FileTest');
gui.mainLayout = uiextras.HBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );

% Create the panels
gui.controlPanel = uiextras.BoxPanel( ...
   'Parent', gui.mainLayout, ...
   'Title', 'MEA Element' );
gui.ViewPanel = uiextras.BoxPanel( ...
   'Parent', gui.mainLayout, ...
   'Title', 'Viewing: ', ...
   'HelpFcn', @onDemoHelp );
% Adjust the main layout
set( gui.mainLayout, 'Sizes', [-1,-3]  );

gui.MeaElement = mysort.mea.H5FileGuiElement(gui.controlPanel);
gui.MeaElement.setFileName(fname);