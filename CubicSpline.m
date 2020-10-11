function [z] = CubicSpline(x,y,v)
%CUBICSPLINE Summary of this function goes here
%   Detailed explanation goes here

% If x or y are column vectors, transpose them into row vectors
[rows,cols] = size(x);
if(rows > 1)
    x = x';
end

[rows,cols] = size(y);
if(rows > 1)
    y = y';
end

% sort x in ascending order and keep the sort index
[x,index] = sort(x);
% rearrange y according to the way x was sorted
y = y(index);

% Code that will evauate each v point
for q = 1:length(v)
    for i = 2:length(x)
        % If this v falls within the ith x range
        if(v(q)<=x(i))
            % Evaluate the (i-1)th function at v
            % z(q) = something
            % break;
        end
    end
end

disp(x);
disp(y);

end

