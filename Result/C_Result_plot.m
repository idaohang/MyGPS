% This function plots the result log from C++ code
% [ tm, ave_p_n, ave_p_e, ave_p_d, shift_n, shift_e, shift_d ]

%load('SHIFT_pos_log20140402_1052.dat');
data = load('SHIFT_pos_log20140407_1334_Max.dat');

true = [0.025, -5.931, 0.051];
tm = data(:,1);
exp_num = length(tm);
pos = data(:,2:4);
shift = data(:,7:9);

% Sats statistic

% Error
one_vec = ones(exp_num,1);
ground_truth = [true(1)*one_vec, true(2)*one_vec, true(3)*one_vec];

pos_error = pos - ground_truth;
bias1 = 0;%mean(error1(:,1));
bias2 = 0;%mean(error1(:,2));
err_norm = zeros(exp_num, 1);
sft_norm = zeros(exp_num, 1);
for i = 1:exp_num
    err_norm(i) = norm([pos_error(i,1)-bias1, pos_error(i,2)-bias2]);
    sft_norm(i) = norm(shift(i,:));
end

ind_sft = find(sft_norm<=2.5);
%err_norm = err_norm(ind_sft);
ind_err = find(err_norm>0.1)
figure(12);
subplot(2,1,1)
[nelements, centers] = hist(err_norm(ind_sft), 30);
bar(centers*100, nelements./sum(nelements)*100, 1);
grid on;
xlabel('Average Horizontal Error (cm)');
ylabel(['Percentage of total ', num2str(sum(nelements)),' trials.'])

subplot(2,1,2)
[nelements, centers] = hist(sft_norm, 30);
bar(centers*1, nelements./sum(nelements)*100, 1);
grid on;
xlabel('NED Shift norm (m)');
ylabel(['Percentage of total ', num2str(sum(nelements)),' trials.'])

% subplot(2,1,2)
% residuals = Mils_Result(:,2:end);
% [II, JJ] = find(abs(residuals)>0.052);
% II
% residuals(~(residuals)) = [];
% residuals = residuals(:);
% %histfit(residuals);hold on;
% [nelements, centers] = hist(residuals, 50);
% bar(centers*1000, nelements./sum(nelements)*100, 1);
% grid on;
% xlabel('Phase measurement residuals (mm)');
% ylabel(['Percentage of total ', num2str(sum(nelements)),' residuals.'])

% max(sft_norm);
% figure(15)
% [nelements, centers] = hist(sft_norm, 30);
% bar(centers*1, nelements./sum(nelements)*100, 1);
% grid on;
% xlabel('NED Shift norm (m)');
% ylabel(['Percentage of total ', num2str(sum(nelements)),' trials.'])
% %else
% %disp('Positioning error is too large!!!')
% %end

max_err = data(:,5);
ind_sft = find(sft_norm<=2.5);
%err_norm = err_norm(ind_sft);
ind_err = find(max_err>0.1)
figure(16)
h(1) = subplot(2,1,1)
[nelements, centers] = hist(max_err, 30);
bar(centers*100, nelements./sum(nelements)*100, 1);
grid on;
%ylabel('Average Horizontal Error norm (cm)')
ylabel('Max Horizontal Error (cm)');

% h(2) = subplot(2,1,2)
% plot(sft_norm(ind_sft))
% grid on;
% ylabel('NED Shift norm (m)')
% xlabel('Index')

% linkaxes(h)
% linkaxes(h,'x')

