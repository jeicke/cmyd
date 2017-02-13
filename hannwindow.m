function w = hannwindow(N)
n = [0:N-1];
w = sin(pi*n/(N-1));