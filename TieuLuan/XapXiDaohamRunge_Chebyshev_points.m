 % Define the Runge function and its derivative
syms x;
f = 1 / (1 + 25*x^2);

% Define range
a = -5;
b = 5;
m = 34;
N = 2*m+1;

% Differentiate the function symbolically
df = diff(f, x);

% Function handles for numerical evaluations
f_func = matlabFunction(f);
df_func = matlabFunction(df);

% Target point for evaluating the derivative
x_target = 0.5;
actual_derivative = df_func(x_target);

% Points for interpolation
% Interpolation points using Chebyshev points
for j = 1:N
    xs(j) = (b+a)/2 + ((b-a)/2)*cos((2*j - 1)*pi/(2*N));
end

ys = f_func(xs);

% Lagrange Interpolating Polynomial
L = 0;
n = length(xs);
for i = 1:n
    Li = 1;
    for j = 1:n
        if i ~= j
            Li = Li * (x - xs(j))/(xs(i) - xs(j));
        end
    end
    L = L + Li * ys(i);
end

% Differentiate the interpolating polynomial
L_prime = diff(L, x);

% Evaluate the derivative of the Lagrange polynomial at the target point
derivative_at_target = double(subs(L_prime, x, x_target));

% Calculate error
error = abs(derivative_at_target - actual_derivative);

fprintf('The actual derivative of the Runge function at x = %.2f is: %.4f\n', x_target, actual_derivative);
fprintf('The approximate derivative from Lagrange interpolation at x = %.2f is: %.4f\n', x_target, derivative_at_target);
fprintf('The error in the derivative approximation is: %.4f\n', error);


