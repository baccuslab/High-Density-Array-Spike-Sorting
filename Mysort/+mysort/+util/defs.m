
function def = defs()
    def.tp = 1;
    def.fn = 2;
    def.fp = 3;
    def.tpo = 4;    % True positive of an overlapping spike
    def.fno = 5;
    def.fpo = 6;
    def.cl  = 7;    % Classification error
    def.clo = 8;    % Classification error of an overlapping spike
    
    def.labelID2String = {'TP', 'FN', 'FP', 'TPO', 'FNO', 'FPO', 'CL', 'CLO'};