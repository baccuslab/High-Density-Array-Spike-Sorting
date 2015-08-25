% function simulateDataTest()
    clear all
    T = [0:10 .5*(0:10); .5*(10:-1:0) 10:-1:0];
    nC = 2;
    Tf = 11;
    
    L = 100;
    noise = randn(nC, L);
    
    gdf = [1 10
           2 15
           1 30
           2 50 
           1 89];
       
    X1 = mysort.util.simulateData(T, gdf, zeros(size(noise)));
    spacer = mysort.plot.mc(X1, 'color', {'g'});
    hold on
    X2 = mysort.util.simulateData(T, gdf, noise);
    mysort.plot.mc(X2, 'figure', 0, 'spacer', spacer, 'color', {'k'});
           