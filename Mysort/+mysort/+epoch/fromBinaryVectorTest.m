% function fromBinaryVectorTest
x = zeros(1,20);

x(1:2) = 1;
x(5:9) = 1;
x(12) = 1;
x(14:16) = 1;
x(18:20) = 1;

epochs = mysort.epoch.fromBinaryVector(x);
% assertEqual(epochs, [1 2
%                      5 9
%                      12 12
%                      14 16
%                      18 20]);
% mysort.plot.figure
% plot(x)
% hold on
% plot(epochs(:,1), 1, 'dr','markersize',15);
% plot(epochs(:,2), 1, '+g','markersize',15);

x = zeros(1,3000);
x(1,608:611) = 1;
epochs = mysort.epoch.fromBinaryVector(x);
% assertEqual(epochs, [608 611]);

