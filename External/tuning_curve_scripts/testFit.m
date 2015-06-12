addpath('C:\Users\frankef\Desktop\SourceCodes\CircStat2011f');
addpath('C:\Users\frankef\Desktop\Visual Coding\tuning_curve_scripts');

vDirRange=0:22.5:359.9;
tuning = [ 0.5 .2 0 1 2 1 0 0.2 0.5 1 3 6 7 6 3 1];
[baseline,maxf1,maxf2,PD,k1,k2,resnorm]=fit_vonmises_fixed_dist(tuning, vDirRange);


fitudir = -2*pi:.1:2*pi;
figure
c1 = mysort.plot.vectorColor(0);
c2 = mysort.plot.vectorColor(1);
thetahat1 = PD*(2*pi/360);
thetahat2 = (PD-180)*(2*pi/360);
[p1 alpha1] = circ_vmpdf(fitudir, thetahat1, k1);
[p2 alpha2] = circ_vmpdf(fitudir, thetahat2, k2);

axes
hold on
plot(vDirRange*(2*pi/360), tuning, 'xk');
plot(alpha1, p1/max(p1)*maxf1, '-', 'color', c1, 'markersize', 12, 'linewidth',2);
plot(alpha2, p2/max(p2)*maxf2, '-', 'color', c2, 'markersize', 12, 'linewidth',2);
