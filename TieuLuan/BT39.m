% Clear variables
clear variables

% Declare symbolic variable x
syms x

% Runge function
f = 1/(1+x^2);

% Choose m: m = 7; m = 10; m = 13
m = 10;

% Interpolation points using Chebyshev points
a = -1;
b = 1;
N = 2*m + 1
% N = 15 => m = 7
% N = 21 => m = 10
% N = 27 => m = 13
for j = 1:N
    X(j) = (b+a)/2 + ((b-a)/2)*cos((2*j - 1)*pi/(2*N))
end

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
N = N(N >= -5 & N <= 5); % Solutions within interpolation interval

% Preallocate GN for efficiency
GN = zeros(length(N), 1);

% Evaluate error at extrema points
for i = 1:length(N)
    GN(i) = subs(G, x, N(i));
end

% No need eval here, due to GN not a symbolic vector
GN = abs(GN);

% Find the maximum error
[Ym, I] = max(GN);
Gmax = subs(abs(G), x, N(I));

% Display Maximum Error
disp('Maximum Absolute Error:');
disp(double(Gmax));