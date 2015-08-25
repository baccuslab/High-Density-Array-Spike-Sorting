
classdef Table < mysort.util.DebuggableClass
    properties
        TH
    end


    methods
        function self = Table(varargin)
            self = self@mysort.util.DebuggableClass(varargin{:});
            P.Units    = 'normalized';
            P.Position = [0 0 1 1];
            P.RearrangeableColumns = 'on';
            [P unresolved] = mysort.util.parseInputs(P, varargin, 'split');
            defP = mysort.util.deflateP(P);
            defUP = mysort.util.deflateP(unresolved);
            self.TH = uitable(defP{:}, defUP{:});   
%             c = {'peter', 1, 'perlich', 2, 2.5
%                  1,  'lala', 'lala', 3, 2};
%             set(self.TH, 'Data', c);
%             set(self.TH, 'ColumnName', {'eins', 'zwei', 'drei', 'vier', 'fünf'});
        end        
    end
end

