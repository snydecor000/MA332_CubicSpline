%% Test the results from the paper: 
% https://ijaers.com/uploads/issue_files/12%20IJAERS-OCT-2017-40-Application%20of%20Cubic%20Spline.pdf

% Data from the paper
x = [7.44 9.30 11.16 13.02 14.88 16.74 18.60]; % Strain of the dogbone in mm
y = [38.76 45.56 49.80 51.89 52.88 52.92 51.68]; % Stress of the dogbone in dkf/mm^2

N = 1000;                               % Number of points the spline will be evaluated at
v = min(x): (max(x)-min(x))/N : max(x); % vector of the x points for the Cubic Spline func

z = CubicSpline(x,y,v);

% Plot it
figure(1);

hold on;
plot(x,y,'o');
plot(v,z);

grid on;
title('Stress-Strain Curve as approximated by a Cubic Spline');
legend('Original Points','Cubic Spline','Location','southeast');
xlabel('Strain    (mm)');
ylabel('Stress    (kgf/mm^2)');
xlim([7 19]);
hold off;