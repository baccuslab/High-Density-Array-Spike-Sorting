M = [ 1  1  1  1  1
      2  3  4  3  2
      14 1  9  4  -1
      -4 2  3  4  2];
  
[i,j,v] = mysort.util.matrixArgMax(M)
assert(i==3 & j==1 & v==14, 'Test failed!');
assert(M(i,j) == v, 'Test failed!');
M = [];
[i,j,v] = mysort.util.matrixArgMax(M)
assert(isempty(i) & isempty(j) & isempty(v), 'Test failed!');

fprintf('Test successful!\n');

M = [-1 0];
[i,j,v] = mysort.util.matrixArgMax(M);
assert(i==1 && j==2)
