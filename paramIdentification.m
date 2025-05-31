%% Скрипт посвящён идентификации параметров ориентационной модели, определении наиболее оптимальных значений в допустимых диапазонах

% загрузка экспериментальных данных, взятых из [W. J. Eldridge, A. Sheinfeld, M. T. Rinehart, and A. Wax, Imaging deformation of adherent cells due to shear stress using quantitative phase imaging, Opt. Lett. 41, 352 (2016).]
% данные содержат 4 зависимости полной деформации от времени выдержки при различных сдвиговых напряжениях
% t_exp - список значений времени выдержки
% eps_exp - список значений полной деформации
% sigm_exp - список сдвиговых напряжений
% eps0_exp - список деформации в начальный момент времени

load('exp-data.mat');

%% Тестовые проверки работы функции CalculateBarrier 
% достигнута непрерывность функции барьера dF при вариации напряжений s0d
% построены профили потенциала F с напряжениями s = -0.6, -0.3, 0.3
close all; clc;
xi = 0.1;

s = [-1:0.01:0.5]';
[dF, ~, ~] = CalculateBarrierF(s, 0.1); 
figure();
plot(s, dF, 'LineWidth', 2);
xlabel('s0d');
ylabel('dF');

s = [-0.6 -0.3 0.3]';
[~, Fset, eRange] = CalculateBarrierF(s,0.1);
figure();
plot(eRange,Fset(1,:),eRange,Fset(2,:),eRange,Fset(3,:),'LineWidth',2);
xlabel('ε0');
ylabel('F-<F>');
legend('s0d = -0.6','s0d = -0.3','s0d = 0.3');

%%
clear all; close all; clc;
load('exp-data.mat');

Gr_range_step = 50;
Ge_range_step = 100e6;
theta_range_step = 10;
L_range_step = 1;

Gr_range = 1:Gr_range_step:1000;
Ge_range = 1e6:Ge_range_step:1e9;
theta_range = (100:theta_range_step:400) * 1.38e-23;
L_range = (1:L_range_step:20) * 1e-6;

Gr_range_length = length(Gr_range);
Ge_range_length = length(Ge_range);
theta_range_length = length(theta_range);
L_range_length = length(L_range);

param_range_defin = false(Gr_range_length, Ge_range_length, theta_range_length, L_range_length);
for i = 1:Gr_range_length
    for j = 1:Ge_range_length
        for k = 1:theta_range_length
            for q = 1:L_range_length
                eps = FindDeformation(t_exp{1}, Gr_range(i), Ge_range(j), theta_range(k), L_range(q), sigm_exp(1), eps0_exp(1));
                param_range_defin(i, j, k, q) = (eps(end) > eps0_exp(1));
            end
        end
    end
end


%%
Gr = 900;
Ge = 1e6;
theta = 300 * 1.38e-23;
L = 10 * 1e-6;
eps = FindDeformation(t_exp{3}, Gr, Ge, theta, L, sigm_exp(3), eps0_exp(3));
plot(t_exp{3}, eps);
median(eps)>eps0_exp

%% Идентификация параметров модели с помощью оптимизации на основе четырёх экспериментальных кривых
close all; clc; clear all;
load('exp-data.mat');

Gr_fixed = 23;                      % Pa, значение эффективного модуля сдвига клетки в представительном объёме
Ge_fixed = 0.4e9;                   % Pa, значение эффективного модуля сдвига сети филаментов в представительном объёме клетки
theta_fixed = 300 * 1.38 * 10^-23;  % J, эффективный температурный фактор
L_fixed = 5e-6;                     % m, длина сегмента филамента
param_fixed_list = [Gr_fixed Ge_fixed theta_fixed L_fixed];

% логический список, указывающий на то, какие параметры будут оптимизированы (1 - оптимизировать, 0 - оставить фиксированное значение ._fixed )
% порядок оптимизируемых параметров соответствует индексам списка [Gr Ge theta L]
which_optim = [1 1 1 1];          

% выполняем множество попыток подбора значений параметров из диапазонов и считаем величину ошибки
N_guesses = 10;
n_param = length(which_optim); % общее число параметров
n_optim_param = sum(which_optim); % число оптимизируемых параметров
guesses = zeros(N_guesses, n_param);
for i = 1:N_guesses
    for j = 1:n_param
        if(which_optim(j) == 0) % если параметр не нужно оптимизировать
            guesses(i,j) = param_fixed_list(j);
        else % если параметр необходимо оптимизировать, то случайным образом выбираем значение параметра из диапазона
            switch j
                case 1 
                    guesses(i,j) = 1e0+(1e3-1e0) * rand;
                case 2
                    guesses(i,j) = 1e6+(1e9-1e6) * rand;
                case 3
                    guesses(i,j) = (100+(400-100)) * 1.38e-23 * rand;
                case 4
                    guesses(i,j) = (0.2+(20-0.2)) * 1e-6 * rand;
                end
        end
    end
end

optim_param_set = zeros(N_guesses, n_param);
rmse_set = zeros(N_guesses, 1);
options = optimset('Display', 'off', 'MaxIter', 1000, 'TolFun', 1e-6, 'TolX', 1e-6);
for i=1:N_guesses
    [optim_param, ~] = fminsearch(@(x) FindObjectiveFunction(x, t_exp, eps_exp, sigm_exp, eps0_exp), ...
        guesses(i,:), options);
    optim_param_set(i,:) = which_optim .* optim_param + ~which_optim .* param_fixed_list;
    rmse_val = FindObjectiveFunction(optim_param_set(i,:), t_exp, eps_exp, sigm_exp, eps0_exp);
    rmse_set(i) = rmse_val;
end

[min_rmse, best_ind] = min(rmse_set);
best_param = optim_param_set(best_ind, :);
%%
p = [71.4459745918108	60963632.2593202	3.65190560044080e-21	6.88115538955539e-06];
[optim_param, ~] = fminsearch(@(x) FindObjectiveFunction(x, t_exp, eps_exp, sigm_exp, eps0_exp), ...
        p, options);
%%
optim_param_set(i,:) = which_optim .* optim_param + ~which_optim .* param_fixed_list;
%% локальные функции
function rmse = FindObjectiveFunction(param_list, t_exp, eps_exp, sigm_exp, eps0_exp) 
Gr = param_list(1); 
Ge = param_list(2);
theta = param_list(3);
L = param_list(4);

error_sum = 0; error_points = 0;
for i=1:length(sigm_exp)
    try
        [~, eps] = ode15s(@(t, eps) deformation_rate(eps, Gr, Ge, theta, L, sigm_exp(i)), ...
            t_exp{i}, eps0_exp(i), odeset('RelTol',1e-8,'AbsTol',1e-10));
        error = (eps_exp{i} - eps).^2;
        error_sum = error_sum + sum(error);
        error_points = error_points + length(error);
    catch
        warning('Error in ode15s for sigma = %.4f Pa', sigm_exp(i));
        error_sum = error_sum + 1e6;
        error_points = error_points + length(t_exp{i});
    end
end
rmse = sqrt(error_sum / error_points);
end

function def = FindDeformation(t, Gr, Ge, theta, L, sigm, eps0) 
[~, def] = ode15s(@(t, eps) deformation_rate(eps, Gr, Ge, theta, L, sigm), ...
    t, eps0, odeset('RelTol',1e-8,'AbsTol',1e-10));
end

function dedt = deformation_rate(eps, Gr, Ge, theta, L, sigm)
% функция скорости деформаций (производной деформаций) для решения ДУ в FindDeformation
xi = 0.1; % xi должно быть = 0.1 (иначе кривые потенциала в диапазоне s0d [-0.6; -0.1; 0.3] не получаются качественно различными)
gam = 16 * 10e-18 * pi * L;
lamd = theta / (xi * gam);
tau0 = 1/(FindOscillationFrequency(L, 3 * Ge));
eps0 = eps - ((sigm - Gr * eps) / ( Ge));
sigm0 = sigm - eps * Gr;
ksiVal = ksiDependence(eps0);
sigm0d = sigm0 ./ (xi * lamd); % напряжения в безразмерном виде
nu = nuCalculation(sigm0d, xi, Ge, tau0);
dedt = (Ge / (nu * (Ge + Gr))) * (sigm - Gr * eps - (theta/(gam)) * ksiVal + lamd * eps0);
end

function f = FindOscillationFrequency(L, E)
% вычисление частоты колебаний сегмента филамента
% L - длина сегмента филамента, E - модуль Юнга филамента
beta = 120.9;
% радиус поперечного сечения равен 4 нм
r = 4e-9; 
roA = 2.53e-14;
I = pi * r^4 / 4;
f = (beta./(2 * pi * L.^2)) .* sqrt(E' .* I / (roA));
end

function ksi = ksiDependence(eps0)
% вычисление величины ksi(eps0)
    ksi = (260 + eps0 .* (123200 + (76250-114200 * eps0) .* eps0))./(12400 + eps0 .* (16770 + eps0 .* (-20670 + (-8820 + eps0) .* eps0)));
end

function [dF, Fset, eRange] = CalculateBarrierF(sod, xi) 
% вычислить величину барьера(-ов) потенциала F, определённом в диапазоне деформаций eRange
% и вывести множество значений потенциала F в диапазоне eRange при различных значениях s0d
% eRange - диапазон деформаций
% s0d - безразмерное значение (список значений) сдвигового напряжения в ориентационном элементе
% xi = theta/(lamda * gamma)
% xStep - шаг по оси деформаций для подсчета данных потенциала
eStep = 0.01;
eRange = -0.49:eStep:0.98;

% выводит множество значений потенциала в диапазоне eRange с набором напряжений s0d
Fset = real(-1.7*log(0.99-eRange) - 114161*log(8822.34-eRange) - 0.9*log(0.5+eRange) - 36.34*log(2.84+eRange)) - (eRange.^2)/(2*xi) - sod.* eRange;
% определяет среднее значение потенциала во всей матрице
mean_val = mean(mean(Fset)); 
% выполнение операции центрирования потенциала F
Fset = Fset - mean_val * ones(1,length(eRange));
dF = zeros(1,length(sod));

% выполнение правил вычисления барьера по заданным профилям F
    for i=1:length(sod)
        F = Fset(i,:);
        max_ind_set = islocalmax(F);
        min_ind_set = islocalmin(F);

        % если кривая профиля F не содержит лок. максимумов и минимумов
        if(numel(eRange(max_ind_set)) + numel(eRange(min_ind_set)) == 0)
            if(F(1) > F(end))
                dF(i) = 0;
            else
                dF(i) = F(end) - F(1);
            end

        % если кривая профиля F содержит лишь один лок. минимум
        elseif(numel(eRange(max_ind_set)) + numel(eRange(min_ind_set)) == 1) 
            min_ind = find(min_ind_set, 1);
            if(eRange(min_ind) > 0)
                dF(i) = 0;
            % барьер определяется разностью потенциалов в 1) точке лок. минимума и 
            % в 2)(!) точке локального минимума первой производной потенциала
            % так достигается непрерывность функции барьера dF при вариации сдвиговых напряжений
            else
                FDeriv = diff(F)./eStep;
                min_deriv_ind = find(islocalmin(FDeriv), 1);
                dF(i) = F(min_deriv_ind) - F(min_ind);
            end

        % если кривая профиля F содержит хотя бы один лок. минимум и максимум
        else
            % определяем индекс с самым левым (по оси eRange) минимумом и максимумом
            max_ind = find(max_ind_set, 1); 
            min_ind = find(min_ind_set, 1); 
            dF(i) = F(max_ind) - F(min_ind);
        end      
    end
end

function nu = nuCalculation(sod, xi, Ge, tau0) 
% вычисление ориентационной вязкости в ОС
% σ0d - безразмерные значения сдвиговых напряжений в ориентационном элементе
% xi = theta/(lamda * gamma)
% Ge - средний модуль сдвига филаментов в ПО клетки
% tau0 - среднее время колебаний филаментов
    nu = Ge * tau0 * exp(CalculateBarrierF(sod, xi));
end