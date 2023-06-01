% simulation settings
tStart = 0;                                     % start time
tEnd = 20;                                      % end time
stepSize = 0.01;                                % size between each element
numEl = length(tStart : stepSize : tEnd);	% number of elements
% the equation
t = linspace(tStart, tEnd, numEl)';	        % create evenly spaced number of elements over the interval [tStart, tEnd] on the x-axis 
y = ((5*exp(-(log(5/3)/5)*(t - 5)) + 10*exp(-(log(2)/5)*t))/2 + ((5*exp(-(log(5/3)/5)*(t - 5)) - 10*exp(-(log(2)/5)*t))/2).*sign(t - 5)).*cos((2*pi/5)*t);
% Create a 2D line plot and specify the line color,  line style, and line width.

fig = figure(1);
fig.Position = [100 100 1000 400];
plot(t, y, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('y');
% title('Attenuated sine wave')

exportgraphics(fig, '..\Figures\Overleaf\theory\exAttenuation.pdf')