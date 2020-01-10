function A = v2a(v,N,M)
%
% A is array of size (N-rows,M-columns)
% v is vector of size (M*N-elements)
%
if (min(size(v)) ~= 1)
    error 'vector is not vector !?!'
end
NM = max(size(v));
if (NM ~= M*N)
    error 'inconsistent dimensions of vector !?!'
end

A = zeros(N,M);
for k=1:NM
    j=mod(k-1,M)+1;
    i=(k-j)/M + 1;
    A(i,j)=v(k);
end