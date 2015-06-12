function diagonalLoadingTest()

    D = mysort.util.diagonalLoading(eye(15), 10);
    assertElementsAlmostEqual(cond(D),1);
    
    D = mysort.util.diagonalLoading(zeros(15), 10);
    assertElementsAlmostEqual(cond(D),1);
    
    r = [1 0.9 0.8];

    CC{1} = toeplitz(r);
    CC{2} = randn(30);
    
    for i=1:length(CC)
        C = CC{i};

        D = mysort.util.diagonalLoading(C, 10);
        assertElementsAlmostEqual(cond(D),10);

        D = mysort.util.diagonalLoading(C, 15.3);
        assertElementsAlmostEqual(cond(D),15.3);    

        D = mysort.util.diagonalLoading(C, 1);
        assertElementsAlmostEqual(cond(D),1);    
    end
    