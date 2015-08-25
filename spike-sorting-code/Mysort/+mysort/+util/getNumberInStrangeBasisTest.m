% mysort.util.getNumberInStrangeBasis(1000, [10 2 4 2 9])
% 
num=mysort.util.getNumberInStrangeBasis(32, [2 2 2 2 2 2])
assert(~any(num~=[1 0 0 0 0 0]), 'Test failed!');

num=mysort.util.getNumberInStrangeBasis(123456789, [10 10 10 10 10 10 10 10 10])
assert(~any(num~=[1 2 3 4 5 6 7 8 9]), 'Test failed!');
