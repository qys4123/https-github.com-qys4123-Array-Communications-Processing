
function Rxx=issp(Rxxm,K)
[M,MM]=size(Rxxm);          %suppose M=5
N=M-K+1;                    %When k=4, N=2, that means the blocks can move two times
J=fliplr(eye(M));           %Inverted diagonal matrix, center symmetry, dimension m * m.
Rxxb=(Rxxm+J*Rxxm.'*J)/2;   %Smoothing the covariance matrix forward and backward (combined with Toeplitz covariance matrix estimation); The signal subspace obtained by conjugate inversion is the same as the original one
Rxx=zeros(K,K);
for i=1:N
Rxx=Rxx+Rxxb(i:i+K-1,i:i+K-1);
end
Rxx=Rxx/N;                  %compute the mean