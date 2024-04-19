% Clear variables
clear variables

% Declare symbolic variable x
syms x

% Runge function
f = 1/(1+x^2);

% Choose m: m = 7; m = 10; m = 13; m = 15
m = 13;

% Interpolation points using linspace in the interval [-5, 5]
%X = linspace(-5, 5, 2*m + 1);

% Interpolation points using linspace in the interval [-2, 2], [-1.5, 1.5]
X = linspace(-1.5, 1.5, 2*m + 1);

% Alternative, you can use for loop here
%for j = 1:2*m+1
%    X(j) = -5+5*(j-1)/m;
%end

% Interpolation values
Y = 1./(1+X.^2);

% Construct the Lagrange interpolating polynomial
PN = 0;
for j = 1:2*m+1
    PN = PN + Y(j) * base_lagrange(X, j, x);
end

% Define the error function and its derivative
G = f - PN;
dG = diff(G, x);

% Solve for zeros of the derivative of G within the interpolation interval
N = double(solve(dG, x));
N = N(imag(N) == 0); % Only consider real solutions
%N = N(N >= -5 & N <= 5); % Solutions within interpolation interval
N = N(N >= -1.5 & N <= 1.5);

% Preallocate GN for efficiency
GN = zeros(length(N), 1);

% Evaluate error at extrema points
for i = 1:length(N)
    GN(i) = subs(G, x, N(i));
end
GN = abs(GN);

% Find the maximum error
[Ym, I] = max(GN);
Gmax = subs(abs(G), x, N(I));

% Display Maximum Error
disp('Maximum Absolute Error:');
disp(double(Gmax));