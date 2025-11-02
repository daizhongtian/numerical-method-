clear; clc;

%% ======== 1. Experiment Parameters ======================================
A      = 1;              
alpha  = 1;              
T      = 2*pi;          
omega0 = 2*pi/T;         
M      = 10;             
N_list = [50,200,1000];  



%% ======== 2. “True” Reference Coefficients (integral) ==================
[a_int, b_int] = fourier_coeffs_builtin(A,alpha,T,M);

%% ======== 3. Compare Rectangle vs. Integral Fourier Coefficients ======
fprintf('\n============== Rectangle vs. Integral Fourier Coefficients ==============\n');
fprintf('  n\t   N=50 (h=%.4f)\tN=200 (h=%.4f)\tN=1000 (h=%.4f)\n',...
        T/N_list(1),T/N_list(2),T/N_list(3));

for n = 0:M
    fprintf('%3d',n);
    for k = 1:numel(N_list)
        [a_rect,b_rect] = fourier_coeffs_rect(A,alpha,T,M,N_list(k));
        if n==0
            val = a_rect(1);
        else
            val = sqrt(a_rect(n+1).^2 + b_rect(n).^2);
        end
        fprintf('\t%12.6e', val);
    end
    fprintf('\n');
end


fprintf('\n-------------- MAE error（rectangle vs. integral） --------------\n');
fprintf('  N\t a_n  MAE\t        b_n  MAE\n');

for k = 1:numel(N_list)                              
    [a_rect, b_rect] = fourier_coeffs_rect( ...     
        A, alpha, T, M, N_list(k));

 
    err_a = mean( abs(a_rect - a_int) );             
    err_b = mean( abs(b_rect - b_int) );             

    fprintf('%5d\t%12.6e\t%12.6e\n', ...          
            N_list(k), err_a, err_b);
end



%% ========== Signal Reconstruction and Visualization ===================
M_plot = [1, 3, 10]; 
N_plot = 200;         

t_plot = linspace(0, T, 1000);
[a_rect_plot, b_rect_plot] = fourier_coeffs_rect(A,alpha,T,max(M_plot),N_plot);
x_orig = x_signal(t_plot, A, alpha, T);

figure('Position',[100 100 1000 600]);
for k = 1:length(M_plot)
    M0 = M_plot(k);
    
    xM_rect = a_rect_plot(1)/2 * ones(size(t_plot));
    xM_int  = a_int(1)/2          * ones(size(t_plot));
    for n = 1:M0
        xM_rect = xM_rect + a_rect_plot(n+1)*cos(n*omega0*t_plot) ...
                           + b_rect_plot(n)  *sin(n*omega0*t_plot);
        xM_int  = xM_int  + a_int(n+1)*cos(n*omega0*t_plot) ...
                           + b_int(n)  *sin(n*omega0*t_plot);
    end
    
    subplot(2,2,k);
    plot(t_plot, x_orig, 'k-', 'LineWidth',1.2); hold on;
    plot(t_plot, xM_rect,'r--','LineWidth',1);
    plot(t_plot, xM_int, 'b-.','LineWidth',1);
    hold off;
    
    title(sprintf('M = %d', M0));
    xlabel('t'); ylabel('x(t)');
    legend('Original signal','Rectangle reconstruction','Integral reconstruction','Location','Best');
    grid on;
end
subplot(2,2,4);
axis off;
text(0,0.5, {
    'Observations and Analysis:'
    '1. For small M (e.g. M=1,3), reconstruction deviates significantly from original;'
    '2. Overshoot/ringing occurs near t=T/2 (Gibbs phenomenon);'
    '3. As M increases (e.g. M=10), reconstruction matches more closely and ringing narrows;'
    }, 'FontSize',12, 'Interpreter','none');

%% extra plot 

M_extra = [15, 20, 25, 30];
N_plot   = 200;
t_plot   = linspace(0, T, 1000);


[a_rect_all, b_rect_all] = fourier_coeffs_rect   (A, alpha, T, max(M_extra), N_plot);
[a_int_all , b_int_all ] = fourier_coeffs_builtin(A, alpha, T, max(M_extra));

for M0 = M_extra

    xM_rect = a_rect_all(1)/2 * ones(size(t_plot));
    xM_int  = a_int_all (1)/2 * ones(size(t_plot));
    for n = 1:M0
        xM_rect = xM_rect ...
                + a_rect_all(n+1).*cos(n*omega0*t_plot) ...
                + b_rect_all(n   ).*sin(n*omega0*t_plot);
        xM_int  = xM_int  ...
                + a_int_all (n+1).*cos(n*omega0*t_plot) ...
                + b_int_all (n   ).*sin(n*omega0*t_plot);
    end
    

    figure('Name',sprintf('Reconstruction M = %d',M0),'NumberTitle','off');
    plot(t_plot, x_orig,   'k-','LineWidth',1.2); hold on;
    plot(t_plot, xM_rect, 'r--','LineWidth',1);
    plot(t_plot, xM_int,  'b-.','LineWidth',1); hold off; grid on;
    title(sprintf('Fourier reconstruction (M = %d)', M0));
    xlabel('t'); ylabel('x(t)');
    legend('Original','Rectangle','Integral','Location','Best');
end

%% ========== Mean-Square Error ε_M Calculation and Plot ================
M_vec = 1:25;
N_err = 200;

eps_rect = zeros(size(M_vec));
eps_int  = zeros(size(M_vec));

Ex2 = integral(@(t) x_signal(t,A,alpha,T).^2, 0, T, ...
               'ArrayValued',true,'RelTol',1e-10,'AbsTol',1e-12);

[a_rect_all,b_rect_all] = fourier_coeffs_rect(A,alpha,T,max(M_vec),N_err);
[a_int_all ,b_int_all ] = fourier_coeffs_builtin(A,alpha,T,max(M_vec));

for idx = 1:numel(M_vec)
    M0 = M_vec(idx);
    Ecap_rect = T*( a_rect_all(1)^2/4 ...
                  + 0.5*sum(a_rect_all(2:M0+1).^2 + b_rect_all(1:M0).^2) );
    eps_rect(idx) = Ex2 - Ecap_rect;

    Ecap_int  = T*( a_int_all(1)^2/4 ...
                  + 0.5*sum(a_int_all(2:M0+1).^2 + b_int_all(1:M0).^2) );
    eps_int(idx)  = Ex2 - Ecap_int;
end

figure;
semilogy(M_vec, eps_rect, 'ro-', 'LineWidth',1.2,'MarkerSize',5); hold on;
semilogy(M_vec, eps_int , 'bs-', 'LineWidth',1.2,'MarkerSize',5);
grid on;
xlabel('Number of harmonics M');
ylabel('Mean-square error \epsilon_M (log scale)');
title('Mean-square error \epsilon_M vs. M');
legend('Rectangle (N=200)','Integral','Location','Best');

N_vals = [50 100 200 400 800 1600 3200 6400];   % ← added four finer meshes
h_vals = T ./ N_vals;                          % step sizes to be plotted

M_coeff = 10;                                  % keep M used for MAE study
errMAE_a = zeros(size(N_vals));
errMAE_b = zeros(size(N_vals));

[a_int_all, b_int_all] = fourier_coeffs_builtin(A, alpha, T, M_coeff);

for k = 1:numel(N_vals)
    [a_rect_k, b_rect_k] = fourier_coeffs_rect(A, alpha, T, M_coeff, N_vals(k));
    errMAE_a(k) = mean(abs(a_rect_k - a_int_all));
    errMAE_b(k) = mean(abs(b_rect_k - b_int_all));
end

% --------- put results in a table and export to CSV ----------------------
maeTable = table(N_vals.', h_vals.', errMAE_a.', errMAE_b.', ...
    'VariableNames', {'N', 'h', 'MAE_a', 'MAE_b'});
disp(maeTable);                         % print nicely in Command Window
writetable(maeTable, 'mae_vs_h.csv');   % creates mae_vs_h.csv in pwd
% -------------------------------------------------------------------------

% ---------- slope fitting (same as before) -------------------------------
pMAE_a = polyfit(log(h_vals), log(errMAE_a), 1);
pMAE_b = polyfit(log(h_vals), log(errMAE_b), 1);

% -------------- refreshed log–log plot -----------------------------------
figure;
loglog(h_vals, errMAE_a, 'ro-', 'LineWidth',1.2,'MarkerSize',6); hold on;
loglog(h_vals, errMAE_b, 'bs-', 'LineWidth',1.2,'MarkerSize',6);

% reference slope‑1 line (optional)
h_ref   = [min(h_vals) max(h_vals)];
refLine = errMAE_a(1)/h_vals(1)^pMAE_a(1) * h_ref.^pMAE_a(1);
loglog(h_ref, refLine, 'k--', 'LineWidth',1);

grid on;
xlabel('Step size h');
ylabel('Mean Absolute Error of Fourier coefficients');
title('MAE vs. step size h (log–log) — extended mesh');
legend(sprintf('a_n  MAE (slope = %.2f)', pMAE_a(1)), ...
       sprintf('b_n  MAE (slope = %.2f)', pMAE_b(1)), ...
       'Reference h^{1}', 'Location','SouthWest');

%% ========================== Local Functions =============================

function y = x_signal(t,A,alpha,T)
% x_signal  periodic signal x(t) over one period (vectorized)
t_mod = mod(t,T);
y     = zeros(size(t_mod));
idx   = (t_mod>0) & (t_mod<=T/2);
y(idx)= A * exp(alpha * t_mod(idx));
end

function I = rect_int(f,a,b,N)
% rect_int  left‑rectangle numerical integration
h = (b-a)/N;
t = a + (0:N-1)*h;
I = h * sum( f(t) );
end

function [a_n, b_n] = fourier_coeffs_rect(A,alpha,T,M,N)
% fourier_coeffs_rect  compute real Fourier coefficients via rectangle rule
omega0 = 2*pi/T;
a_n    = zeros(1,M+1);
b_n    = zeros(1,M);

% a0 term
I0     = rect_int(@(t)x_signal(t,A,alpha,T), 0, T, N);
a_n(1) = 2/T * I0;

% a_n, b_n (n=1..M)
for n = 1:M
    fcos = @(t) x_signal(t,A,alpha,T) .* cos(n*omega0*t);
    fsin = @(t) x_signal(t,A,alpha,T) .* sin(n*omega0*t);
    a_n(n+1) = 2/T * rect_int(fcos, 0, T, N);
    b_n(n)   = 2/T * rect_int(fsin, 0, T, N);
end
end

function [a_n, b_n] = fourier_coeffs_builtin(A,alpha,T,M)
% fourier_coeffs_builtin  compute real Fourier coefficients via Matlab integral
omega0 = 2*pi/T;
a_n    = zeros(1,M+1);
b_n    = zeros(1,M);
opts   = {'ArrayValued',true,'RelTol',1e-10,'AbsTol',1e-12};

% a0 term
a_n(1) = 2/T * integral(@(t)x_signal(t,A,alpha,T), 0, T, opts{:});

% a_n, b_n (n=1..M)
for n = 1:M
    a_n(n+1) = 2/T * integral(@(t)x_signal(t,A,alpha,T).*cos(n*omega0*t), ...
                              0, T, opts{:});
    b_n(n)   = 2/T * integral(@(t)x_signal(t,A,alpha,T).*sin(n*omega0*t), ...
                              0, T, opts{:});
end
end
