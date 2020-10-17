function [z] = CubicSpline(x,y,v)
%CUBICSPLINE Summary of this function goes here
%   Detailed explanation goes here
%
%   Cubic Spline Eq: f0(x) = a0 + b0(x-x0) + c0(x-x0)^2 + d0(x-x0)^3
%                             where x is within [x0 and x1]
%
%
%   Format of the coefficient matrix for the cublic spline equations:
%
%   [ a0 b0 c0 d0 ;
%     a1 b1 c1 d1 ;
%     a2 b2 c2 d2 ;
%     ............;
%     an bn cn dn ]
%

% If x or y are column vectors, transpose them into row vectors
[rows,cols] = size(x);
if(rows > 1)
    x = x';
end

[rows,cols] = size(y);
if(rows > 1)
    y = y';
end

% Sort x in ascending order and keep the sort index
[x,index] = sort(x);
% Rearrange y according to the way x was sorted
y = y(index);

% Initialize the output vector and the A matrix
z = zeros(1,length(v));
A = zeros(length(x)-2,length(x)-2);

% Calculate the distance between the x points
deltax = diff(x);

% Calculate the diagonal of the A matrix
for i = 1 : length(deltax)-1
    mainDiagA(i) = 2*(deltax(i) + deltax(i+1));
end

% Construct the A matrix 
A(1,1) = mainDiagA(1);
for i = 2 : length(deltax)-1
    A(i,i) = mainDiagA(i);
    A(i,i-1) = deltax(i);
    A(i-1,i) = deltax(i);
end

B = zeros(1,length(y)-2);
for i = 2 : length(y)-1
    B(i-1) = 6*((y(i+1)-y(i))/deltax(i)-(y(i)-y(i-1))/deltax(i-1));
end

% Solve the system of equations for the second derivative values
s = TDMS(A,B);
% Set the second derivative to be 0 for the first and last points
s = [0 s 0];
% Loop through and solve for the coeffs for each cubic equation based on
% the second derivative values
coeffs = zeros(length(x)-2,4);
for i = 1 : length(y)-1
    coeffs(i,1) = y(i);
    coeffs(i,2) = (y(i+1)-y(i))/deltax(i) -...
                  (deltax(i)*(s(i+1)+2*s(i)))/6;
    coeffs(i,3) = s(i)/2;
    coeffs(i,4) = (s(i+1)-s(i))/(6*deltax(i));
end

% Code that will evauate each v point
for q = 1:length(v)
    for i = 2:length(x)
        % If this v falls within the ith x range
        if(v(q)<=x(i))
            dx = v(q) - x(i-1);
            % Evaluate the (i-1)th function at v
            z(q) = coeffs(i-1,1) + ...
                   coeffs(i-1,2)*dx +...
                   coeffs(i-1,3)*dx^2 +...
                   coeffs(i-1,4)*dx^3;
            break;
        end
    end
end
figure(1)
hold on;
plot(x,y,'o');
plot(v,z);
hold off;
end

