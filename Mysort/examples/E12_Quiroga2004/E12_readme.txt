This is the most complex example from the data analysis and result handling
point of view. Here, we use the "mysort.util.AnalysisHandler" class to 
manage the whole analysis process. This class is initialised with the
following parameters:
    1) a name
    2) an evaluation function
    3) a list (cell array) of parameter names
    4) a list (cell array) of parameter sets for every parameter name
    5) a path where the results will be stored.

The analysis handler will then create all possible parameter-value
combinations. Once the analysis process is started, the evaluation
function 2) will be called for all those parameter combinations and the
results stored. 

Another analysis handler can now be initialised and instead of providing
a save path, we can pass an analysis handler to it. This will allow the 
second instance to build its parameter combinations against all those of 
the first class.

Different evaluations can now be derived from the first class without the 
need to recompute the evaluations of the first class.

In this example a third analysis handler is than using the second to 
do the final evaluation.

