function z = CubicSpline(x,y,v)
%CUBICSPLINE takes points on the xy plane, and synthesizes a natural cubic
%spline which is then evaluated at a number of x coordinates to produce a
%vector of output values
%   
%   z = CUBICSPLINE(x,y,v)
%
%   x and y are parallel row or column vectors that define the (x,y) points
%   on which the cubic spline is created. The (x,y) points can be in any
%   particular order
%
%   v is an optional argument of x coordinates where the cublic spline will
%   be evaluated. It too can be a row or column vector
%
%   z is the resulting vector from evaluating the cublic spline at the 
%   locations in v
%
%   Natural Cubic Splines are created by defining cubic piecewise 
%   equations from a list of (x,y) points. In total, there will be n-1
%   total cubic equations where n is the number of (x,y) points. The
%   following criteria are used in order to constrain the cubic equations
%   and create a natural style cubic spline:
%
%   1.  Each piecewise cubic equation must pass through the beginning and 
%       end points of the piecewise segment 
%   2.  The first derivative of any two adjacent cubic equations should be
%       equal
%   3.  The second derivative of any two adjacent cubic equations should be
%       equal.  The second derivative at the end points should be 0.
%
%   This set of criteria fully constrains a system of equations which can
%   be solved to find the coefficients of each cubic equation in the
%   piecewise set. 
%
%   Piecewise Cubic Eq: f(x) = a + b(x-x0) + c(x-x0)^2 + d(x-x0)^3
%
%                               ( where x is within [x0 and x1] )
%
%   Fortunately, we are able to rearrange and define this system as a 
%   Tridiagonal Matrix System which lets us use some very efficient methods
%   to solve for the coefficients. TDMS() is used to solve this system
%
%   Each Piecewise Cubic Eq. within the Cublic Spline has 4 coefficients,
%   and these values are stored in a matrix to make evaluating the cubic 
%   spline simpler:
%
%   [ a0 b0 c0 d0 ;
%     a1 b1 c1 d1 ;
%     a2 b2 c2 d2 ;
%     ............;
%     an bn cn dn ]
%
%   Eq. 0 applies to all points [x0,x1)
%   Eq. 1 applies to all points [x1,x2)
%   etc...
%
%   2020 Cory Snyder

% If x or y are column vectors, transpose them into row vectors
[rows,cols] = size(x);
if(rows > 1)
    x = x';
end
[rows,cols] = size(y);
if(rows > 1)
    y = y';
end

% If the v argument doesn't exist, then default it to the values in x
if (~exist('v','var'))
    v = x;
end

% If v is a column vector, transpose it into a row vector
[rows,cols] = size(v);
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

% Calculate the main diagonal of the A matrix
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

% Preare the b vector to pass to TDMS 
b = zeros(1,length(y)-2);
for i = 2 : length(y)-1
    b(i-1) = 6*((y(i+1)-y(i))/deltax(i)-(y(i)-y(i-1))/deltax(i-1));
end

% Solve [A][x] = [b] where x represents the second derivative values
s = TDMS(A,b);

% Set the second derivative to be 0 for the first and last points
s = [0 s 0];

% Loop through and solve for the coeffs for each cubic equation based on
% the second derivative values
coeffs = zeros(length(x)-2,4);
for i = 1 : length(y)-1
    coeffs(i,1) = y(i);
    coeffs(i,2) = (y(i+1)-y(i))/deltax(i) - (deltax(i)*(s(i+1)+2*s(i)))/6;
    coeffs(i,3) = s(i)/2;
    coeffs(i,4) = (s(i+1)-s(i))/(6*deltax(i));
end

% Evauates the spline at each v point
for q = 1:length(v)
    % Loop through the values of x to find the correct cubic equation
    for i = 2:length(x)
        % If this v falls within the this (ith) x range
        if(v(q)<=x(i))
            dx = v(q) - x(i-1);
            % Evaluate the (i-1)th cubic equation at v
            z(q) = coeffs(i-1,1) + ...
                   coeffs(i-1,2)*dx +...
                   coeffs(i-1,3)*dx^2 +...
                   coeffs(i-1,4)*dx^3;
            break;
        end
    end
end

% Plot the results 
close 1;
figure(1)
hold on;
title('Resulting Cubic Spline');
plot(x,y,'o');
plot(v,z);
legend('Original Points','Cubic Spline');
xlabel('x');
ylabel('y    ','Rotation',0);
hold off;
end

