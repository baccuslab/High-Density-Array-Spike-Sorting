%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Example without subsample shift in X (see below for other)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subTf = 21;
nC = 2;
T = [[0 0 0 0:5 4:-1:-5 -4:0]/5 0 0 0 0 0 0 [0:5 4:-1:-5 -4:0]/5 0 0 0
      0 0 0 sin((0:subTf-1)*2*pi/subTf) 0 0 0 0 0 0 sin((0:subTf-1)*2*pi/subTf) 0 0 0];
Tf = size(T,2)/nC;
offset = 4;
subT = T(:, mysort.wf.vSubIdx(Tf, nC, offset:offset+subTf-1));
X = repmat(T,100,1);
IDs = repmat([1:2]', 100,1);
tau = repmat([0; 0; repmat([-3:3]', 14, 1)], 2, 1);
% tau = zeros(200,1);
X = mysort.wf.vShift(X, nC, tau, 1);
nX = X +randn(size(X))*.02;
%%
figure;
subplot(5,1,1)
plot(T');
title('original full templates with zero paddings');
subplot(5,1,2)
plot(subT');
title('original templates without zero paddings');
subplot(5,1,3)
plot(X');
title('instances of templates');
subplot(5,1,4)
plot(nX');
title('noisy instances of templates');

%%


[D maxTaus TupShiftDown FupShiftDown tauRange] = mysort.wf.vTemplateMatching(nX, subT, nC, offset, 'maxShift', 5, 'upsample', 5);
tauRange
EF = diag(T*T');                     % compute energies
DISCR = D - .5 * repmat(EF', size(X,1), 1)  + log(.01);  % compute botm dicriminant

[D_ IDs_] = max(DISCR,[],2);
tauMax_ = maxTaus(:,1);
tauMax_(IDs_==2) = maxTaus(IDs_==2,2);

[IDs IDs_ tauMax_ tau]

%%
figure;
plot(DISCR(1:2:end,1), DISCR(1:2:end,2), 'b.');
hold on
plot(DISCR(2:2:end,1), DISCR(2:2:end,2), 'r.');
plot([0 6], [0 6], 'k-');

%%
offset = 4;
subT = T(:, mysort.wf.vSubIdx(Tf, nC, offset:offset+subTf-1));
[D maxTaus TupShiftDown FupShiftDown tauRange] = mysort.wf.vTemplateMatching(nX, subT, nC, offset, 'maxShift', 0, 'upsample', 5);
EF = diag(T*T');                     % compute energies
DISCR = D - .5 * repmat(EF', size(X,1), 1)  + log(.01);  % compute botm dicriminant

[D_ IDs_] = max(DISCR,[],2);
tauMax_ = maxTaus(:,1);
tauMax_(IDs_==2) = maxTaus(IDs_==2,2);

%[IDs IDs_ tauMax_ tau]

plot(DISCR(1:2:end,1), DISCR(1:2:end,2), 'g.');
hold on
plot(DISCR(2:2:end,1), DISCR(2:2:end,2), 'm.');
plot([0 6], [0 6], 'k-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Example WITH subsample shift in X 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subTf = 21;
nC = 2;
T = [[0 0 0 0 0 0:5 4:-1:-5 -4:0]/5 0 0 0 0 0 0 0 0 0 0 [0:5 4:-1:-5 -4:0]/5 0 0 0 0 0
      0 0 0 0 0 sin((0:subTf-1)*2*pi/subTf) 0 0 0 0 0 0 0 0 0 0 sin((0:subTf-1)*2*pi/subTf) 0 0 0 0 0];
Tf = size(T,2)/nC;
offset = 6;
subT = T(:, mysort.wf.vSubIdx(Tf, nC, offset:offset+subTf-1));
X = repmat(T,100,1);
IDs = repmat([1:2]', 100,1);
tau = -4 + 8.*rand(200,1);
% tau = zeros(200,1);
X = mysort.util.shiftRowsInterpolated(X, tau, nC);
nX = X +randn(size(X))*.01;
%%
figure;
subplot(5,1,1)
plot(T');
title('original full templates with zero paddings');
subplot(5,1,2)
plot(subT');
title('original templates without zero paddings');
subplot(5,1,3)
plot(X');
title('instances of templates');
subplot(5,1,4)
plot(nX');
title('noisy instances of templates');

%%

[D maxTaus TupShiftDown FupShiftDown tauRange] = mysort.wf.vTemplateMatching(nX, subT, nC, offset, 'maxShift', 5, 'upsample', 5);
tauRange
EF = diag(T*T');                     % compute energies
DISCR = D - .5 * repmat(EF', size(X,1), 1)  + log(.01);  % compute botm dicriminant

[D_ IDs_] = max(DISCR,[],2);
tauMax_ = maxTaus(:,1);
tauMax_(IDs_==2) = maxTaus(IDs_==2,2);

[IDs IDs_ tauMax_ tau]

%%
figure;
ah = subplot(1,2,1)
plot(DISCR(1:2:end,1), DISCR(1:2:end,2), 'b.');
hold on
plot(DISCR(2:2:end,1), DISCR(2:2:end,2), 'r.');
plot([0 6], [0 6], 'k-');

subplot(1,2,2);
plot(tau(1:2:end), tauMax_(1:2:end), 'bx');
hold on
plot(tau(2:2:end), tauMax_(2:2:end), 'rx');
plot([-5 5], [-5 5], 'k-');
xlabel('real tau');
ylabel('estimated tau');

%%
offset = 4;
subT = T(:, mysort.wf.vSubIdx(Tf, nC, offset:offset+subTf-1));
[D maxTaus TupShiftDown FupShiftDown tauRange] = mysort.wf.vTemplateMatching(nX, subT, nC, offset, 'maxShift', 0, 'upsample', 1);
EF = diag(T*T');                     % compute energies
DISCR = D - .5 * repmat(EF', size(X,1), 1)  + log(.01);  % compute botm dicriminant

[D_ IDs_] = max(DISCR,[],2);
tauMax_ = maxTaus(:,1);
tauMax_(IDs_==2) = maxTaus(IDs_==2,2);

%[IDs IDs_ tauMax_ tau]
axes(ah)
plot(DISCR(1:2:end,1), DISCR(1:2:end,2), 'g.');
hold on
plot(DISCR(2:2:end,1), DISCR(2:2:end,2), 'm.');
plot([0 6], [0 6], 'k-');


