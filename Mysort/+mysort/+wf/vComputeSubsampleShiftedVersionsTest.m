
T = [0 1 2 3 0 0 1 2 3 0
     0 3 2 1 0 0 3 2 1 0];
nC = 2;

    
vT = mysort.wf.vComputeSubsampleShiftedVersions(T, nC, 3);

figure;
subplot(2,1,1);
plot(vT(:,:,1)');
subplot(2,1,2);
hold on
for i=1:size(vT,3)
    plot(vT(:,:,i)', 'color', mysort.plot.vectorColor(i))
end

%%
xrange = 0:.4:2*pi;
T = [sin(xrange) sin(xrange)
     cos(xrange*2) cos(xrange*2)];
nC = 2;
    
vT = mysort.wf.vComputeSubsampleShiftedVersions(T, nC, 20);

figure;
subplot(2,1,1);
plot(vT(:,:,1)');
subplot(2,1,2);
hold on
for i=1:size(vT,3)
    plot(vT(:,:,i)', 'color', mysort.plot.vectorColor(i))
end
