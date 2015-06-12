
configList(1).name = 'TestConfig1';
xcords = repmat(1:10, 10, 1);
configList(1).ELPOS = [xcords(:) repmat((1:10)', 10, 1) (1:100)'];
configList(1).selectedElectrodes = 1:10;

configList(2).name = 'TestConfig2';
xcords = repmat(1:100, 100, 1);
configList(2).ELPOS = [xcords(:) repmat((1:100)', 100, 1) (1:10000)'];
configList(2).selectedElectrodes = 100:110;

activeConfig = 1;
save('configList.mat', 'configList', 'activeConfig', '-v7.3');