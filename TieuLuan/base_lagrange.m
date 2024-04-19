% define base lagrange function
function f=base_lagrange(X,k,x)
N = length(X);
f = 1;
for i=1:N
    if i ~= k
        f = f*(x-X(i))/(X(k)-X(i));
    end
end