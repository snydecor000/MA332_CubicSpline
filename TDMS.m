function x = TDMS(A,b)
%TDMS uses the special LU factorization and special backsolving properties 
%of Tridiagonal Matrix Systems to solve Ax = b
%
%   x = TDMS(A,b)
%
%   A is a square Tridiagonal Matrix and b is a row or column vector with 
%   the same length/width as A.  
%
%   First, this function breaks A down into L and U. Because A is a TDM, we
%   only had to eliminate the leading entry in row i+1 when doing
%   elimination on column i. This significantly simplifies the LU
%   factoriztion algorithm.
%
%   With A broken into L and U, we now need to solve LUx = b. This can be
%   broken into two steps:
%
%   1.  LUx = L(Ux) = Ly = b
%   2.  Ux = y
%
%   Both of these steps become incredibly easy when the unique shape of L
%   and U is exploited. Because A was a TDM, L and U will have the form:
%
%   [  1    0    0    0  ;
%     l21   1    0    0  ;
%      0   l32   1    0  ;
%      0    0   l43   1  ]
%
%   This allows for a fast forward substitution algorithm to be used
%
%   2020 Cherise McMahon

% Initialize U and L and find the size of A
[n,z] = size(A);
U = A;
L = eye(n);

% Use the special LU factorization algorithm because A is a TDM
for i = 1 : n-1
    L(i+1,i) = U(i+1,i)/U(i,i);
    U(i+1,:) = U(i+1,:) - (L(i+1,i)*U(i,:));
end

% Solve Ly=b
y = zeros(1,n);
y(1) = b(1);
for i = 2 : n
    y(i) = b(i) - (L(i,i-1)*y(i-1));
end

% Solve Ux=y
x = zeros(1,n);
x(n)=y(n)/U(n,n);
for i = n-1:-1:1
   x(i) = (y(i) - U(i,i+1)*x(i+1))/U(i,i);
end
