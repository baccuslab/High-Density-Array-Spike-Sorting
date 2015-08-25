
classdef ProcessTimer < handle
    
    properties
        nItems
        nDigits
        tElapsedThis
        tElapsedSmooth
        tLastTime
        tStart
        alpha
        i
    end
    
    methods
        %%% ----------------CONSTRUCTOR---------------------------     
        function self = ProcessTimer(nItems, alpha)
            if nargin < 2
                alpha = .3;
            end
            self.alpha = alpha;
            self.nItems = nItems;
            self.nDigits = floor(log10(nItems))+1;
            self.tElapsedSmooth = [];
            self.tElapsedThis = [];
            self.tStart = [];
            self.i = 0;
        end
        
        %%% ------------------------------------------------------
        function next(self, i)
            if nargin == 1
                self.i = self.i+1;
            else
                self.i = i;
            end
            if isempty(self.tStart)
                self.tStart = tic;
            else
                self.tElapsedThis = toc(self.tLastTime);
            end
            if isempty(self.tElapsedSmooth) || self.tElapsedSmooth < .1
                self.tElapsedSmooth = self.tElapsedThis;
            else
                self.tElapsedSmooth = (1-self.alpha)*self.tElapsedSmooth + ...
                                      self.alpha*self.tElapsedThis;
            end
            self.tLastTime = tic;
        end
        
        %%% ------------------------------------------------------
        function showProgress(self)
            disp(self.getProgressString());
        end
        
        %%% ------------------------------------------------------
        function str = getProgressString(self)
            i = self.i;
            N = self.nItems;
            t = self.tElapsedSmooth;
            p = round(100*(i-1)/N);
            str1 = ['Processing %' num2str(self.nDigits) 'd/%' num2str(self.nDigits) 'd (%3d%% done'];
            if isempty(t)
                str = sprintf([str1 ')'],i,N,p);
                return
            end
            r = N-(i-1);
            T = r*t;
            minu = floor(T/60);
            sek  = ceil(mod(T,60));
            str = sprintf([str1 ', est: % 4d min:%02d sec remaining):'],i,N,p,minu,sek);                    
        end
        %%% ------------------------------------------------------
        function showTotal(self)
            disp('Total: ')
            disp(toc(self.tStart));
        end
    end
end