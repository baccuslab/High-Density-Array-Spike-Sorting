% From http://code.google.com/p/scemua/source/browse/trunk/tools/ggobiexport.m?r=2
function ggobiexport(M,varargin)
%
% <a href="matlab:web(fullfile(scemroot,'html','ggobiexport.html'),'-helpbrowser')">View HTML documentation for this function in the help browser</a>    
%
% % Usage:
% %   ggobiexport(M,...)
% % 
% % Authorized options:
% %   xmlFile
% %   dataSetName
% %   description
% %   varNames
% %   glyphStr
% %   color
% %   formatStr

xmlFile = 'tmp.xml';
dataSetName='dataset name';
description = 'no description available';
nVars = size(M,2);
for k=1:nVars
    varNames{k,1}=['var',num2str(k,'%02d')];
end

glyphStr = 'fr 1';
color = 5;
formatStr = repmat('%10.6f',[1,nVars]);

parsepairs(varargin)


xmlStr=['<?xml version="1.0"?>',char(10),...
'<!DOCTYPE ggobidata SYSTEM "ggobi.dtd">',char(10),...
'',char(10),...
'<ggobidata count="1">',char(10),...
'<data name="',dataSetName,'">',char(10),...
'<description>',char(10),...
description,char(10),...
'</description>',char(10),...
'<variables count="',num2str(size(M,2)),'">',char(10)];

for k=1:nVars
    xmlStr = [xmlStr,...
    '<variable name="',varNames{k},'"/>',char(10)];
end

xmlStr=[xmlStr,...
'</variables>',char(10),...
'<records count="',num2str(size(M,1)),'" glyph="',glyphStr,...
'" color="',num2str(color),'">',char(10)];

fid = fopen(xmlFile,'wt');
fprintf(fid, '%s', xmlStr);
fclose(fid);

xmlStr='';

fid = fopen(xmlFile,'at+');

disp([mfilename,':  0 %'])
for k=1:size(M,1)

    if ismember(k,ceil(size(M,1)*[5:5:95]/100))
        disp([sprintf('\b\b\b\b\b'),sprintf('%2d',100*k/size(M,1)),' %'])
        fprintf(fid, '%s', xmlStr);        
        xmlStr='';
    end
    xmlStr=[xmlStr,'<record label="',num2str(k),'">',char(10),...
            sprintf(formatStr,M(k,:)),char(10),'</record>',char(10)];
end
disp(sprintf(repmat('\b',[1,numel(mfilename)+8])))

xmlStr=[xmlStr,'</records>',char(10),...
'</data>',char(10),...
'</ggobidata>'];


fprintf(fid, '%s', xmlStr);
fclose(fid);

% % % % % % % % % % % % % % % % % % % % % % % % % % 
function parsepairs(C)

AuthorizedOptions={'xmlFile',...
                   'dataSetName',...
                   'description',...
                   'varNames',...
                   'glyphStr',...
                   'color',...
                   'formatStr'};

for k = 1:2:length(C(:))
    if ~strcmp(C{k}, AuthorizedOptions)
        s=dbstack;
        error(['Unauthorized parameter name ' 39 C{k} 39 ' in ' 10,...
            'parameter/value passed to ' 39 s(end-1).name 39 '.']);
    end
    evalin('caller',[C{k},'=varargin{',num2str(k+1),'};'])
end