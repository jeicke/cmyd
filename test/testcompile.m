A = GPUsingle(randn(5));
B = GPUsingle(randn(5));
% A and B are dummy variables
GPUcompileStart('myfun','-f', A, B)
R1 = exp(A);
R2 = floor(B);
GPUcompileStop(R1,R2)