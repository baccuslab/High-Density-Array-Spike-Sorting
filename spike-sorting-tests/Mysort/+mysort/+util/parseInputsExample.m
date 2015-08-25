    P = [];
    P.bla1 = 'bla';
    P.bla2 = 'bla';
    
    P2.bla2 = 'Peter';
    P2.bla3 = 'hallo';
   
    [P_ uP] = mysort.util.parseInputs(P, {P2, 'aber', 'hallo'});  
    P_
    uP
    
    [P_ uP] = mysort.util.parseInputs(P, {P2, 'aber', 'hallo'}, 'merge');
    P_
    uP
    
%     fprintf('THIS SHOULD THROW AN ERROR:\n');
%     [P_ uP] = mysort.util.parseInputs(P, 'test', {P2, 'aber', 'hallo'}, 'error');
