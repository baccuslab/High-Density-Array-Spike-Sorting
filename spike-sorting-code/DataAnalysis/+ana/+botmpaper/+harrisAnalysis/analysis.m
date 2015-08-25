classdef analysis < ana.analysisWrapperInterface
    properties
        
    end
    methods(Abstract, Access=public)
        [figureHandles figureNames] = makeFigures_(self, name, result, P);
    end    
    methods(Abstract, Access=protected)
        result = runTrial_(self, name);
    end
    methods
        %------------------------------------------------------------------
        function self = analysis(analysisName)
            [path outpath] = ana.harris.datapath();
 
            dataInputPath = path;
            dataOutPutPath = outpath;
            figurePath = ana.botmpaper.harrisAnalysis.figureOutPath();
            self = self@ana.analysisWrapperInterface(analysisName, dataInputPath, dataOutPutPath, figurePath);  
        end
        %------------------------------------------------------------------
        function validNames = getValidTrialNames(self)
            info = ana.harris.info();
            validNames = info(:,2);
        end
    end
end
