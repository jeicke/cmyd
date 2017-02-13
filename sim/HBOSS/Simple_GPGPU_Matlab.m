%% Matlab(R) General-Purpose computation on Graphics Processing Units (GPGPU) 
% 
% The following code requires GPUmat from www.gp-you.org
% Copyright 2009 gp-you.ch.
%
% Follow the instructions (on top of the command window): press "next"
% to continue the example or "stop" to stop it. 

%%
% First we have to start the GPU world. You should install GPUmat from
% gp-you.org in order to run this example. 
%
% The GPU device is initialized by typing GPUstart. The GPUmat installation
% folder should be included in the path to run this command. If you receive
% an error please check the GPUmat User Guide for help.  

%% Start the GPU 
GPUstart

%%
% GPUs have their own memory. To perform operations on the GPU, first we
% have to create variables on the GPU memory. This is done in MATLAB(R) by
% using the GPUsingle class. A GPUsingle is a single precision floating
% point variable. There are several ways to create a GPUsingle, the easiest
% way is to build one from a MATLAB(R) array. 
%
% So, first we create a MATLAB(R) array Ah and then we create a GPU
% variable A from it.

%% Creation of a variable on GPU memory
% Ah is a MATLAB(R) array. Ah is on HOST memory. HOST is your PC, where the
% GPU is installed. Ah is created using the single() command, because we
% want to use single precision floating point variables.
% 
% A is a GPU variable created from Ah. A is on GPU memory. A is a GPU
% single precision floating variable.

Ah = single(rand(5,5)) % Ah is on HOST memory
A  = GPUsingle(Ah)     % A  is on GPU memory

%% 
% Operations on GPU variables are executed on the GPU. Please check the
% GPUmat manual for available functions.

%% Execution on the GPU
% We create a second GPU variable B and we add it to A, storing the result
% in C, which will be also a GPU variable. The sum is performed on the GPU.

B = GPUsingle(rand(5,5)); % B is on GPU memory
C = A+B                   % executed on GPU. C is also on GPU memory

%% 
% To perform the same operation on the CPU we create a second MATLAB(R)
% array Bh and add it to Ah, storing the result in Ch. Bh is created from
% the GPU variable B using the casting function single(B), which converts
% the GPUsingle B into the MATLAB(R) array Bh. 
%
% Please note that the command used to execute the code on the GPU (A+B) is
% the same as the command used to execute the code on the CPU (Ah+Bh). The
% difference is just the type of variables.

Bh = single(B);  % Bh is on CPU memory, created from B using single()
Ch = Ah+Bh       % executed on CPU

%%
% Again...
% Ah and Bh are on HOST memory -> Ah+Bh is executed on the CPU
% A and B are on GPU memory    -> A+B is executed on the GPU

C = A+B     % executed on GPU. C is also on GPU memory
Ch = Ah+Bh  % executed on CPU

%% Some more complicated operation

C  = exp(A).*B + B.^2    % executed on GPU
Ch = exp(Ah).*Bh + Bh.^2 % executed on CPU

%% Convert GPU variables into Matlab variables
% We use the single() function to convert a GPUsingle into a MATLAB(R)
% array

Ah = single(A); % A is on GPU memory, Ah is on HOST memory
Bh = single(B); % B is on GPU memory, Bh is on HOST memory
Ch = single(C); % C is on GPU memory, Ch is on HOST memory

%% More ways to create variables on the GPU...
% Create a vector using the colon function. It is equivalent to using the
% MATLAB(R) colon function.

G = colon(1,10,50,GPUsingle)   % G is on GPU memory
Gh = single(1:10:50)           % Gh is on HOST memory

%% Using the colon function...
H = colon(50,-10,1,GPUsingle)  % H is on GPU memory
Hh = single(50:-10:1)          % Hh is on HOST memory

%% Using the colon function...
I = colon(1,0.1,1.5,GPUsingle) % I is on GPU memory
Ih = single(1:0.1:1.5)         % Ih is on HOST memory

%% Complex numbers 
% Using G and H to create a sequence of complex numbers on the GPU
L = sqrt(-1) * G          % L is on GPU memory, complex
M = sqrt(-1) * pi * 2 * H % M is on GPU memory, complex

%% Using vertical concatenation to create GPU variables
N = [zeros(5,1,GPUsingle);colon(1,1,5,GPUsingle)';zeros(5,1,GPUsingle)] 

%% Selecting some elements from GPU variable N using a vector of indexes
idx = GPUsingle([6 7 8]); % idx is the index vector on GPU memory
O = N(idx)

%% FFT transform 
% The following code will show how to create two GPU variables A and B and
% to perform 1D and 2D FFT on them.

%% Create two arrays A and B on the GPU 
% Initialize them with random numbers
A = GPUsingle(rand(1,100));   % A is on GPU memory
B = GPUsingle(rand(100,100)); % B is on GPU memory

%% 1D FFT (GPU)
% The 1D FFT is executed on the GPU, because A is a GPU variable
FFT_A = fft(A); % executed on GPU

%% FFT, GPU Vs.CPU execution
% A is a GPUsingle           -> fft(A)  is executed on GPU
% Ah is a MATLAB(R) variable -> fft(Ah) is executed on CPU

Ah = single(A);   % A is on GPU memory, Ah is on HOST memory
FFT_A  = fft(A);  % executed on GPU
FFT_Ah = fft(Ah); % executed on CPU

%% 2D FFT (GPU)
FFT_B = fft2(B); % executed on GPU

%% FFT2, GPU Vs.CPU execution
% B is a GPUsingle           -> fft2(B)  is executed on GPU
% Bh is a MATLAB(R) variable -> fft2(Bh) is executed on CPU

Bh = single(B);    % B is on GPU memory, Bh is on HOST memory
FFT_B  = fft2(B);  % executed on GPU
FFT_Bh = fft2(Bh); % executed on CPU


%% GPU Vs. CPU speed up
% Now we create to big arrays and multiply them element-by-element. The CPU
% execution time is compared to the GPU execution time, and the speed up is
% calculated
%
% GPUsync blocks MATLAB(R) execution until all the GPU tasks have been
% completed. It is required only for benchmarks.
%
% NOTE: The GPU requires some warm up. It might be that you get better
% speed up the second time you run this example. 

A = GPUsingle(rand(4000,4000)); % A is on GPU memory
B = GPUsingle(rand(4000,4000)); % B is on GPU memory

%% GPU benchmark

% Warm up GPU
C = A.*B;

% Run test
tic;
C = A.*B;
GPUsync; 
% GPUsync required only on benchmarks to make sure
% that GPU operations have finished
timegpu = toc;

%% CPU benchmark
% First convert GPU variables into Matlab variables Ah and Bh
Ah = single(A);
Bh = single(B);
tic;
Ch = Ah.*Bh;
timecpu = toc;

%% calculate speed up
speedup = timecpu/timegpu

%% Complex numbers
% The following example shows the previous benchmark with complex numbers.

A = colon(1,1,4e6,GPUsingle)*sqrt(-1);  % A is on GPU memory, complex
B = colon(4e6,-1,1,GPUsingle)*sqrt(-1); % B is on GPU memory, complex

%% GPU benchmark
tic;
C = A.*B;
GPUsync; 
% GPUsync required only on benchmarks to make sure
% that GPU operations have finished
timegpu = toc;

%% CPU benchmark
% First convert GPU variables into Matlab variables Ah and Bh
Ah = single(A);
Bh = single(B);
tic;
Ch = Ah.*Bh;
timecpu = toc;

%% calculate speed up
speedup = timecpu/timegpu







