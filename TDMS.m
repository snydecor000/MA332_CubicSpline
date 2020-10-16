function [x] = TDMS(A,b)
%TDMS Summary of this function goes here
%   Detailed explanation goes here

% Initialize U and L and find the size of A
[n,z] = size(A);
U = A;
L = eye(n);

% Use special LU factorization because A is a TDM
for i = 1 : n-1
    L(i+1,i) = U(i+1,i)/U(i,i);
    U(i+1,:) = U(i+1,:) - L(i+1,i)*U(i,:);
end

% Initialize the output x
x = zeros(1,n);

% Solve the L matrix for the output values using forward substitution
x(1) = b(1);
for i = 2 : n
    x(i) = b(i) - L(i,i-1)*x(i-1);
end
