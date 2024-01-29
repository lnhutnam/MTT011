clear all;
syms x y z;

w = [x;y;z];

f = 2*x^2 - x + y^2 - z;
g = 32*x^2 - y^2 -20*z;
h = y^2 - 14*x*z;

H = [f;g;h];

J = [diff(f,x) diff(f,y) diff(f,z); diff(g,x) diff(g,y) diff(g,z); diff(h,x) diff(h,y) diff(h,z)];

w1 = subs(w,{x,y,z},{0.5, 1.0, 0.0});

J1 = subs(J,{x,y,z},{0.5, 1.0, 0.0});

w2 = w1 - inv(J1) * subs(H,{x,y,z},{0.5, 1.0, 0.0})


J2 = subs(J,{x;y;z},w2)

