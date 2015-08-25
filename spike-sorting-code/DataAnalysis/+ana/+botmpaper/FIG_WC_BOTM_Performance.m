% This matrix contains evaluation results
%    WC           BOTM     BOTM+SIC
% det cla tot |det cla tot |   det cla tot   |   det cla tot
%  1   2   3    4   5    6      7   8   9      10  11  12
X = [...
921	1	922     76	2	78      13	4	17     11   2   13
236	5	241     45	3	48      5	2	7      4    2   6
374	5	379     59	5	64      5	0	5      8    0   8
999	12	1011	61	9	70      9	3	12     9    3   12
174	3	177     51	19	70      2	1	3      2    2   4
193	10	203     39	21	60      7	1	8      6    2   8
184	45	229     43	36	79      4	2	6      4    6   10
637	306	943     85	93	178     6	2	8      6    7   13
274	0	274     37	17	54      3	17	20      2   18  20
201	41	242     39	50	89      18	11	29      18  12  30
217	81	298     46	132	178     11	17	28      9   17  26
405	651	1056	45	278	323     19	10	29      20  10  30
183	1	184     45	34	79      9	8	17      8   9   17
157	8	165     33	38	71      4	7	11      5   9   14
193	443	636     51	153	204     5	8	13      8   9   17
492	1462 1954	104	386	490     7	14	21      36  1124    1160
5840 3074 8914	859	1276 2135	127	107	234     156 1232 1388];


D = load('QuirogaBenchmarkDetails.mat');

%%
n_spikes_per_file = sum(D.nSpikesPerNeuron,2);

detIdx = [1 4 7 10];
claIdx = [2 5 8 11];

Neval = length(detIdx);

TP_det = repmat(n_spikes_per_file, [1 Neval])-X(1:end-1, detIdx);
TP_cla = repmat(n_spikes_per_file, [1 Neval])-X(1:end-1, claIdx);

P_det = round(1000 * TP_det./repmat(n_spikes_per_file, [1 Neval]))/10;
P_cla = round(1000 * TP_cla./repmat(n_spikes_per_file, [1 Neval]))/10;
P_tot = round(10*(P_det+P_cla)/2)/10;

det_total = round(1000 * sum(TP_det)./(ones(1, Neval)*sum(n_spikes_per_file)))/10
cla_total = round(1000 * sum(TP_cla)./(ones(1, Neval)*sum(n_spikes_per_file)))/10
(det_total + cla_total)/2

%%
fname = 'QuirogaBigTable.txt';
try
    delete(fname);
end
XX = [];
for i=1:4
    NDet = X(1:end-1,(i-1)*3+1);
    PDet = 100*(n_spikes_per_file-NDet)./n_spikes_per_file ;
    NCL  = X(1:end-1,(i-1)*3+2);
    PCL  = 100*(n_spikes_per_file-NCL)./n_spikes_per_file;
    NTot = NDet + NCL;
    assert(~any(NTot ~= X(1:end-1,(i-1)*3+3)), 'e');
    PTot = (PDet+PCL)/2;
    
    XX = [XX NDet PDet NCL PCL NTot PTot];
end
XX = [XX; sum(XX,1)];
XX(end, 2:2:end) = 100* (sum(n_spikes_per_file)-XX(end, 1:2:end)) ./ sum(n_spikes_per_file);
for i=1:17
    str = [sprintf('%d (%.1f), ', XX(i,:)) '\n'];
    mysort.util.appendToFile(fname, str);
end


%%
figure
bar(P_det)