function parseInputsTest()

    P2.bla2 = 'Peter';
    P2.bla3 = 'hallo';

    Pneu = lala('otherP', P2, 'bla1', 'kalle');
    
    assertEqual(Pneu.bla1, 'kalle')
    assertEqual(Pneu.bla2, 'Peter')
    f = @() Pneu.bla3;
    assertExceptionThrown(f, 'MATLAB:nonExistentField');
    
    function P = lala(varargin)
        P.bla1 = 'bla';
        P.bla2 = 'bla';
        
        P = mysort.util.parseInputs(P, 'test', varargin);        
    end

end