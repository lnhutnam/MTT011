function bfd = backwardFiniteDifference(f, x, h, m, n)
    % Hàm tính gia số lùi hữu hạn (backward finite difference)
    % Inputs:
    %   f: function handle
    %   x: điểm đánh giá
    %   h: kích thước bước nhảy
    %   m: bậc đạo hàm (1 với đạo hàm cấp một, 2 với đạo hàm cấp hai)
    %   n: bậc sai số (1 for O(h), 2 for O(h^2))
    % Output:
    %   bfd: xấp xỉ gia số lùi hữu hạn của đạo hàm

    % Kiểm tra đầu vào
    if m > 2 || n > 2
        error('m and n must be <= 2');
    end

    coefficients = [
        NaN, NaN,  1, -1;      % Đạo hàm cấp một, O(h)
        NaN, 3/2, -2, 1/2;     % Đạo hàm cấp hai, O(h^2)
    ];
    coefficients(:, :, 2) = [
        NaN,  1, -2,  1;       % Đạo hàm cấp một, O(h^2)
         2, -5,  4, -1;        % Đạo hàm cấp hai, O(h^2)
    ];

    % Hệ số cho bậc đạo hàm và bậc sai số
    c = coefficients(n, :, m);
    c = c(~isnan(c));

    % Tính toán gia số lùi hữu hạn
    bfd = 0;
    for i = 1:length(c)
        bfd = bfd + c(i) * f(x - (i-1) * h);
    end
    bfd = bfd / (h^m);
end

function ffd = forwardFiniteDifference(f, x, h, m, n)
    % Hàm tính gia số tiến hữu hạn (forward finite difference)
    % Tương tự như Hàm tính gia số lùi hữu hạn (backward finite difference)
    % Tham số đầu vào:
    %   f: function handle
    %   x: điểm đánh giá
    %   h: kích thước bước nhảy
    %   m: bậc đạo hàm (1 với đạo hàm cấp một, 2 với đạo hàm cấp hai)
    %   n: bậc sai số (1 for O(h), 2 for O(h^2))
    % Kết quả đầu ra:
    %   ffd: xấp xỉ gia số tiến hữu hạn của đạo hàm

    if m > 2 || n > 2
        error('m and n must be <= 2');
    end

    coefficients = [
        -1,  1, NaN, NaN;      % Đạo hàm cấp một, O(h)
        -3/2, 2, -1/2, NaN;    % Đạo hàm cấp hai, O(h^2)
    ];
    coefficients(:, :, 2) = [
         1, -2,  1, NaN;       % Đạo hàm cấp một, O(h^2)
         2, -5,  4, -1;        % Đạo hàm cấp hai, O(h^2)
    ];

    c = coefficients(n, :, m);
    c = c(~isnan(c));

    ffd = 0;
    for i = 1:length(c)
        ffd = ffd + c(i) * f(x + (i-1) * h);
    end
    ffd = ffd / (h^m);
end

function fcd = firstCenteredDerivative(f, x, h)
    % Tính toán đạo hàm cấp một trung tâm (first centered derivative)
    fcd = (f(x+h) - f(x-h)) / (2*h);
end

function scd = secondCenteredDerivative(f, x, h)
    % Tính toán đạo hàm cấp hai trung tâm (second centered derivative)
    scd = (f(x+h) - 2*f(x) + f(x-h)) / (h^2);
end

% Định nghĩa biến và hàm symbolic
syms x;
f(x) = sin(x^2);

% Chuyển đổi hàm symbolic thành hàm MATLAB để tính toán số học
func1 = matlabFunction(f);

% Initialize parameters for step size variation
% Khởi tạo tham số cho kích thước bước nhảy
pow = 0.6; % yếu tổ mở rộng (expansion)
lim = 50;  % cấp số cao nhất (limitation)
i = 1:lim;
arr = @(x) double(power(pow,x));
hs = arr(i); % sinh ra mảng bước nhảy hs = [h^1, h^2, ..., h^lim]

% Đỉnh đánh giá đạo hàm
x_eval = 2;

% Tính giá trị đúng của đạo hàm bậc nhất và đạo hàm bậc hai tại x
df_actual   = double(subs(diff(f, 1), x, x_eval)); % Đạo hàm cấp một
d2f_actual  = double(subs(diff(f, 2), x, x_eval)); % Đạo hàm cấp hai

% Khởi tạo mảng để tìm lỗi và xấp xỉ
errors = zeros(10, length(hs));
approximations = zeros(1, 10);

% Bảng lưu trữ lỗi tối thiểu cho mỗi phương pháp
minErrors = zeros(10, 2);

% Lặp lại từng giá trị h để tính toán sai số cho các phương pháp lấy sai phân khác nhau
for j = 1:length(hs)
    h = hs(j); % Bước nhảy hiện tại

    % Tính gần đúng cho đạo hàm bậc nhất và bậc hai
    % Centered first derivative - O(h^2)
    approximations(1)  = firstCenteredDerivative(func1, x_eval, h);

    % Centered second derivative - O(h^2)
    approximations(2)  = secondCenteredDerivative(func1, x_eval, h);

    % Forward Finite first derivative - O(h)
    approximations(3)  = forwardFiniteDifference(func1, x_eval, h, 1, 1);

    % Forward Finite first derivative - O(h^2)
    approximations(4)  = forwardFiniteDifference(func1, x_eval, h, 1, 2);

    % Forward Finite second derivative - O(h)
    approximations(5)  = forwardFiniteDifference(func1, x_eval, h, 2, 1);

    % Forward Finite second derivative - O(h^2)
    approximations(6)  = forwardFiniteDifference(func1, x_eval, h, 2, 2);

    % Backward Finite first derivative - O(h)
    approximations(7)  = backwardFiniteDifference(func1, x_eval, h, 1, 1);

    % Backward Finite first derivative - O(h^2)
    approximations(8)  = backwardFiniteDifference(func1, x_eval, h, 1, 2);

    % Backward Finite second derivative - O(h)
    approximations(9)  = backwardFiniteDifference(func1, x_eval, h, 2, 1);

    % Backward Finite second derivative - O(h^2)
    approximations(10) = backwardFiniteDifference(func1, x_eval, h, 2, 2);

    % Tính toán và lưu trữ các sai số tương đối cho mỗi phép tính gần đúng
    for i = 1:10
        if mod(i, 2) == 1 % Chỉ số lẻ cho đạo hàm bậc nhất
            errors(i, j) = abs(approximations(i) / df_actual - 1);
        else % Chỉ số chẵn cho đạo hàm bậc hai
            errors(i, j) = abs(approximations(i) / d2f_actual - 1);
        end
    end
end

% Vẽ đồ thị sai số cho phép tính gần đúng đạo hàm trung tâm bậc nhất
figure;
loglog(hs, errors(1,:), 'r-*');
hold on;
[minError, idxMinError] = min(errors(1,:));
loglog(hs(idxMinError), errors(1,idxMinError), "ro", 'MarkerSize', 7, 'LineWidth', 3);
set(gca, 'XDir', 'reverse');
xlabel("h (step size)");
ylabel("Relative Error (%)");
title("Error vs. step size h for central finite difference of first derivative", 'FontSize', 12);
legend("f^{(1)}(x) + O(h^2)", "Min error");
minErrors(1, :) = [hs(idxMinError), errors(1, idxMinError)];
hold off;
%-------------------------------------------------------------------
% Vẽ đồ thị sai số cho phép tính gần đúng đạo hàm trung tâm bậc hai
figure;
loglog(hs,errors(2,:),'r-*')
hold on;
[~,idx] = min(errors(2,:));
loglog(hs(idx),errors(2,idx),"r o",'MarkerSize',7,'LineWidth',3)
set(gca,'XDir','reverse')
xlabel("h (step size)")
ylabel("Relative Error (%)")
[t,~] = title("Error vs length of h for central finite difference of second derivative",' ');
t.FontSize = 12;
legend("f^{(2)}(x) +O(h^2)","min error")
tabErrors(2,1)= hs(idx);
tabErrors(2,2)= errors(2,idx);
hold off
%-------------------------------------------------------------------
% Vẽ đồ thị sai số cho phép tính gần đúng sai phân tiến cho đạo hàm bậc một
figure;
loglog(hs,errors(3:4,:),"-*")
hold on
[~,idx] = min(errors(3,:));
loglog(hs(idx),errors(3,idx),"b o",'MarkerSize',7,'LineWidth',3)
tabErrors(3,1)= hs(idx);
tabErrors(3,2)= errors(3,idx);


[~,idx] = min(errors(4,:));
loglog(hs(idx),errors(4,idx),"r o",'MarkerSize',7,'LineWidth',3)
tabErrors(4,1)= hs(idx);
tabErrors(4,2)= errors(4,idx);
set(gca,'XDir','reverse')
xlabel("h (step size)")
ylabel("Relative Error (%)")
[t,~] = title("Error vs length of h for Forward Finite Difference of first derivative",' ');
t.FontSize = 12;
legend("f^{(1)}(x) +O(h)","f^{(1)}(x) +O(h^2)","min error (h)","min error (h^2)")
hold off
%-------------------------------------------------------------------
% Vẽ đồ thị sai số cho phép tính gần đúng sai phân tiến cho đạo hàm bậc hai
figure;
loglog(hs,errors(5:6,:),"-*")
hold on
[~,idx] = min(errors(5,:));
loglog(hs(idx),errors(5,idx),"b o",'MarkerSize',7,'LineWidth',3)
tabErrors(5,1)= hs(idx);
tabErrors(5,2)= errors(5,idx);

[~,idx] = min(errors(6,:));
loglog(hs(idx),errors(6,idx),"r o",'MarkerSize',7,'LineWidth',3)
set(gca,'XDir','reverse')
xlabel("h (step size)")
ylabel("Relative Error (%)")
[t,s] = title("Error vs length of h for Forward Finite Difference of second derivative", ' ');
t.FontSize = 12;
legend("f^{(2)}(x) +O(h)","f^{(2)}(x) +O(h^2)","min error (h)","min error (h^2)")
tabErrors(6,1)= hs(idx);
tabErrors(6,2)= errors(6,idx);
hold off
%-------------------------------------------------------------------
% Vẽ đồ thị sai số cho phép tính gần đúng sai phân lùi cho đạo hàm bậc một
figure;
loglog(hs,errors(7:8,:),"-*")
hold on
[~,idx] = min(errors(7,:));
loglog(hs(idx),errors(7,idx),"b o",'MarkerSize',7,'LineWidth',3)
tabErrors(7,1)= hs(idx);
tabErrors(7,2)= errors(7,idx);

[~,idx] = min(errors(8,:));
loglog(hs(idx),errors(8,idx),"r o",'MarkerSize',7,'LineWidth',3)
set(gca,'XDir','reverse')
xlabel("h (step size)")
ylabel("Relative Error (%)")
title("Error vs length of h for Backward Finite Difference of first derivative",' ')
legend("f^{(1)}(x) +O(h)","f^{(1)}(x) +O(h^2)","min error (h)","min error (h^2)")
tabErrors(8,1)= hs(idx);
tabErrors(8,2)= errors(8,idx);
hold off
%-------------------------------------------------------------------
% Vẽ đồ thị sai số cho phép tính gần đúng sai phân lùi cho đạo hàm bậc hai
figure;
loglog(hs,errors(9:10,:),"-*")
hold on
[~,idx] = min(errors(9,:));
loglog(hs(idx),errors(9,idx),"b o",'MarkerSize',7,'LineWidth',3)
tabErrors(9,1)= hs(idx);
tabErrors(9,2)= errors(9,idx);

[~,idx] = min(errors(10,:));
loglog(hs(idx),errors(10,idx),"r o",'MarkerSize',7,'LineWidth',3)
set(gca,'XDir','reverse')
xlabel("h (step size)")
ylabel("Relative Error (%)")
title("Error vs length of h for Backward Finite Difference of second derivative",' ')
legend("f^{(2)}(x) +O(h)","f^{(2)}(x) +O(h^2)","min error (h)","min error (h^2)")
hold off
%-------------------------------------------------------------------
pre = 2:2:14; % Bậc chính xác: 2, 4, 6, 8, 10, 12, 14
m = 1; % Bậc đạo hàm
hs = arr(1:30);

% Khởi tạo ma trận độ lỗi
error = zeros(length(pre), length(hs));

for idx = 1:length(pre)
    n = pre(idx);
    n_coefs = 2*floor((m+1)/2)-1+n; p = (n_coefs-1)/2;
    A = power(-p:p, (0:2*p)');
    b = zeros(2*p+1,1);
    b(m+1) = factorial(m);
    c = A\b; % Các hệ số cho xấp xỉ sai phân hữu hạn

    for j = 1:length(hs)
        h = hs(j);
        k = -p;
        ffd_val = 0;
        for cof = c'
            ffd_val = ffd_val + cof * func1(x_eval + k * h);
            k = k + 1;
        end
        ffd_val = ffd_val / (h^m);

        % Tính toán và lưu trữ sai số tương đối cho độ chính xác và kích thước bước hiện tại
        error(idx, j) = abs((ffd_val - df_actual) / df_actual);
    end
end

figure
l(1) = loglog(hs,error(1,:),"b-*");
hold on
l(2) = loglog(hs,error(2,:),"r-*");
l(3) = loglog(hs,error(3,:),"g-*");
l(4) = loglog(hs,error(4,:),"m-*");
l(5) = loglog(hs,error(5,:),"c-*");
l(6) = loglog(hs,error(6,:),"y-*");
l(7) = loglog(hs,error(7,:),"k-*");
set(gca,'XDir','reverse');
[~,idx] = min(error(1,:));
loglog(hs(idx),error(1,idx),"b o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(2,:));
loglog(hs(idx),error(2,idx),"r o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(3,:));
loglog(hs(idx),error(3,idx),"g o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(4,:));
loglog(hs(idx),error(4,idx),"m o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(5,:));
loglog(hs(idx),error(5,idx),"c o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(6,:));
loglog(hs(idx),error(6,idx),"y o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(7,:));
loglog(hs(idx),error(7,idx),"k o",'MarkerSize',8,'LineWidth',3)
xlabel("h (step size)")
ylabel("Relative Error (%)")
title("\fontsize{10}Error of central finite difference for the first derivative for differents precisions",' ')
legend(l(1: 7),'h^2','h^4','h^6','h^8','h^{10}','h^{12}','h^{14}');
hold off
%-------------------------------------------------------------------
pre = 2:2:14; % Bậc chính xác: 2, 4, 6, 8, 10, 12, 14
m = 2; % Bậc đạo hàm
hs = arr(1:30);

% Khởi tạo ma trận độ lỗi
error = zeros(length(pre), length(hs));

for idx = 1:length(pre)
    n = pre(idx);
    n_coefs = 2*floor((m+1)/2)-1+n; p = (n_coefs-1)/2;
    A = power(-p:p, (0:2*p)');
    b = zeros(2*p+1,1);
    b(m+1) = factorial(m);
    c = A\b; % Các hệ số cho xấp xỉ sai phân hữu hạn

    for j = 1:length(hs)
        h = hs(j);
        k = -p;
        ffd_val = 0;
        for cof = c'
            ffd_val = ffd_val + cof * func1(x_eval + k * h);
            k = k + 1;
        end
        ffd_val = ffd_val / (h^m);

        % Tính toán và lưu trữ sai số tương đối cho độ chính xác và kích thước bước hiện tại
        error(idx, j) = abs((ffd_val - d2f_actual) / d2f_actual);
    end
end

figure
l(1) = loglog(hs,error(1,:),"b-*");
hold on
l(2) = loglog(hs,error(2,:),"r-*");
l(3) = loglog(hs,error(3,:),"g-*");
l(4) = loglog(hs,error(4,:),"m-*");
l(5) = loglog(hs,error(5,:),"c-*");
l(6) = loglog(hs,error(6,:),"y-*");
l(7) = loglog(hs,error(7,:),"k-*");
set(gca,'XDir','reverse');
[~,idx] = min(error(1,:));
loglog(hs(idx),error(1,idx),"b o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(2,:));
loglog(hs(idx),error(2,idx),"r o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(3,:));
loglog(hs(idx),error(3,idx),"g o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(4,:));
loglog(hs(idx),error(4,idx),"m o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(5,:));
loglog(hs(idx),error(5,idx),"c o",'MarkerSize',8,'LineWidth',3)
[~,idx] = min(error(6,:));
loglog(hs(idx),error(6,idx),"y o",'MarkerSize',8,'LineWidth',3)
[minV,idx] = min(error(7,:));
loglog(hs(idx),error(7,idx),"k o",'MarkerSize',8,'LineWidth',3)
xlabel("h (step size)")
ylabel("Relative Error (%)")
title("\fontsize{10}Error of central finite difference for the second derivative for differents precisions",' ')
legend(l(1: 7),'h^2','h^4','h^6','h^8','h^{10}','h^{12}','h^{14}');
hold off
