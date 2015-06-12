function shiftRowsTest()

    X = [0 0 1 0 0
         0 0 1 0 0];
    Y1 = [0 0 1 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 1 0 0];
    Y1_ = mysort.util.shiftRows(X, [0 5], 0);
%     Y1_ == Y1
   assertEqual(Y1,Y1_);
    
    Y2 = [0 0 1 0 0
          0 0 0 0 0];      
    Y2_ = mysort.util.shiftRows(X, [0 5], 1);
%     Y2_ == Y2
   assertEqual(Y2,Y2_);

    Y3 = [0 0 0 0 0 0 0 1 0 0
          0 0 1 0 0 0 0 0 0 0];
    Y3_ = mysort.util.shiftRows(X, [0 -5], 0);
%     Y3_ == Y3
   assertEqual(Y3,Y3_);      

    Y4 = [0 0 1 0 0 
          0 0 0 0 0 ];
    Y4_ = mysort.util.shiftRows(X, [0 -5], 1);
%     Y4_ == Y4
   assertEqual(Y3,Y3_);  


    