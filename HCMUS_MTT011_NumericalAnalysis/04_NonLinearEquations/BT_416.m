% Tren len chung ta da giai quyet


% clear all;
% syms x y z;

% w = [x;y;z];

% f = 2*x^2 - x + y^2 - z;
% g = 32*x^2 - y^2 -20*z;
% h = y^2 - 14*x*z;

% H = [f;g;h];

% J = [diff(f,x) diff(f,y) diff(f,z); diff(g,x) diff(g,y) diff(g,z); diff(h,x) diff(h,y) diff(h,z)];

% w1 = subs(w,{x,y,z},{0.5, 1.0, 0.0});

% J1 = subs(J,{x,y,z},{0.5, 1.0, 0.0});

% w2 = w1 - inv(J1) * subs(H,{x,y,z},{0.5, 1.0, 0.0})

% J2 = subs(J,{x;y;z},w2)


function BT416()
    % Khai bao cac bien
    syms x y z;

    % Dinh nghia he phuong trinh phi tuyen
    f = 2*x^2 - x + y^2 - z;
    g = 32*x^2 - y^2 -20*z;
    h = y^2 - 14*x*z;

    % Dao ham rieng
    % Doi voi ham f
    dfdx = diff(f, x);
    dfdy = diff(f, y);
    dfdz = diff(f, z);

    % Doi voi ham g
    dgdx = diff(g, x);
    dgdy = diff(g, y);
    dgdz = diff(g, z);

    % Doi voi ham h
    dhdx = diff(h, x);
    dhdy = diff(h, y);
    dhdz = diff(h, z);

    % Diem khoi tao
    x0 = 0.5;
    y0 = 1.0;
    z0 = 1.0;

    % Tolerance, va so luot lap lon nhat
    tolerance = 1e-10;
    max_iterations = 25;

    % Bat dau qua trinh lap
    for n = 1:max_iterations % Matlab bat dau tai 1, bruh
        % Danh gia gia tri ham va dao ham tai [x_0, y_0, z_0]
        f_value = subs(f, [x, y, z], [x0, y0, z0]);
        g_value = subs(g, [x, y, z], [x0, y0, z0]);
        h_value = subs(h, [x, y, z], [x0, y0, z0]);

        % Tinh toan ma tran Jacobian
        J = [subs(dfdx, [x, y, z], [x0, y0, z0]), subs(dfdy, [x, y, z], [x0, y0, z0]), subs(dfdz, [x, y, z], [x0, y0, z0]); ...
             subs(dgdx, [x, y, z], [x0, y0, z0]), subs(dgdy, [x, y, z], [x0, y0, z0]), subs(dgdz, [x, y, z], [x0, y0, z0]); ...
             subs(dhdx, [x, y, z], [x0, y0, z0]), subs(dhdy, [x, y, z], [x0, y0, z0]), subs(dhdz, [x, y, z], [x0, y0, z0])]; % Jacobian matrix

        % Tinh toan sai so delta x, delta y, delta z
        delta = J \ [-f_value; -g_value; -h_value]; 

        % Cap nhat
        x0 = x0 + delta(1);
        y0 = y0 + delta(2);
        z0 = z0 + delta(3);

        % Kiem tra dieu kien hoi tu
        if norm([delta(1); delta(2); delta(3)], inf) < tolerance
            break;
        end
    end

    % Tra ra ket qua
    fprintf('Approximate root: x = %.10f, y = %.10f, z = %.10f\n', x0, y0, z0);
    fprintf('Number of iterations: %d\n', n);

end