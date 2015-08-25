ffile = '/links/groups/hima/recordings/collaborations/Antonia/NHP150416/p1c1/Trace_id1274_2015-04-16T14_57_45_0.stream.ntk';

siz = 10000;
ntk=initialize_ntkstruct(ffile, 'nofilters');
[high_density_data ntk] = ntk_load(ntk, siz);
nC = size(ntk.data, 1);

chunkSize = 500000;
FRN_Numbers = ntk.frame_no;
while ~ntk.eof
    [high_density_data ntk] = ntk_load(ntk, chunkSize);
    FRN_Numbers = [FRN_Numbers ntk.frame_no];
end
FRN_Numbers = double(FRN_Numbers);
%%
figure
plot(1:length(FRN_Numbers), FRN_Numbers -( (1:length(FRN_Numbers))+FRN_Numbers(1)-1))

length(FRN_Numbers)/20000