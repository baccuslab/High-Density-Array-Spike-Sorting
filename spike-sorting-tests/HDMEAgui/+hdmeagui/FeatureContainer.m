classdef FeatureContainer < handle
    properties
       % internal members
       featureList
       featureNames
    end
    
    methods
        %------------------------------------------------------------------
        function self = FeatureContainer(featureList)
            self.featureList = featureList;
        end
    end
end