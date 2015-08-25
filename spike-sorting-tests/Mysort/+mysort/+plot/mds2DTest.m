%% Create Electrode Array
nEls = 25;
range = linspace(0, 1, sqrt(nEls));
[EX, EY] = meshgrid(range, range);
Epos = [EX(:) EY(:)];

%% Create Positions of Templates on 2D Array as regular Grid
nNeurons = 4*4;
range = linspace(0, 1, sqrt(nNeurons));
[NX, NY] = meshgrid(range, range);
Npos = [NX(:) NY(:)];
Npos = Npos + randn(size(Npos))*.05;
IDs = 1:nNeurons;
C =  mysort.plot.vectorColor(IDs);

%% Create Templates
motherwf = [1:9 8:-1:1];
motherwf = motherwf / norm(motherwf);
L = length(motherwf);
T = zeros(L, nEls, nNeurons);
A = zeros(nNeurons, nEls);
for n=1:nNeurons
    pn = Npos(n,:);
    for i=1:nEls
        pe = Epos(i,:);
        d = norm(pn-pe);
        A(n,i) = normpdf(d, 0, .15);
        T(:,i,n) = motherwf * A(n,i);
    end
end

%% Compute distances between templates
vT = mysort.wf.t2v(T);
Deuk = pdist(vT);
Dint = normcdf(Deuk, 0, 1.5);

figure
plot(Deuk(:), Dint(:), '.')


%%
FAC = 200;
fh = mysort.plot.figure([800 600]);
mysort.plot.templates2D(T, Epos*FAC, 5, [], 'fh', fh, 'IDs', IDs);
hold on
plot(FAC*Epos(:,1), FAC*Epos(:,2), 'k.', 'markersize', 12);
scatter(FAC*Npos(:,1), FAC*Npos(:,2), [], C, 'linewidth', 2);
set(gca, 'xlim', [-FAC/3 FAC+FAC/3], 'ylim', [-FAC/3 FAC+FAC/3]);

%%
opt = statset('Display', 'iter');
Dint_ = Dint;
%Dint_(Dint_>.93) = nan;
Dint_ = Dint_ * 1000;
stress = inf;
for i=1:20
    [Y_,stress_,disparities_] = mdscale(Dint_, 2, 'start', 'random', 'Options', opt);
    if stress_ < stress
        Y = Y_;
        stress = stress_;
        disparities = disparities_;
    end
end
[Y_,stress_,disparities_] = mdscale(Dint_, 2, 'start', Npos, 'Options', opt);
%%
figure
subplot(1,2,1)
scatter(Y_(:,1), Y_(:,2), [], C, 'linewidth', 2)

subplot(1,2,2)
scatter(Y(:,1), Y(:,2), [], C, 'linewidth', 2)