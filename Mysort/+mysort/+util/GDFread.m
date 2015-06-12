
function X = GDFread(fname, delimiter)
    if nargin < 2
        delimiter = ' ';
    end
    fileID=fopen(fname);
    gdfdata=textscan(fileID,['%d' delimiter '%d']);
    [units]=gdfdata{1};
    [samples]=gdfdata{2};
    fclose(fileID);
    X = double([units samples]);
    
%     fileID=fopen(path);
%     gdfdata=textscan(fileID,'%3d%2d %d');
%     [tetrodenr]=gdfdata{1};
%     [neuronID]=gdfdata{2};
%     [Sampletime]=gdfdata{3};
%     fclose(fileID);
%     index=[1:length(Sampletime)];