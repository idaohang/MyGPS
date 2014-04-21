% This function plot the results after revision 73
% [start_time, start_time+buff_size, ...
%    n, sqrt(norm(N1)^2/length(N1)), error1, std1, sqrt(norm(N2)^2/length(N2)), ...
%    error2, std2, sqrt(norm(N3)^2/length(N2)), error3, std3 ];

exp_num = length(Result(:,1));

%% Sats statistics
sat_num = Result(:,3);
ind_r_10 = find(sat_num==10);
ind_r_9 = find(sat_num==9);
ind_r_8 = find(sat_num==8);
ind_r_7 = find(sat_num==7);
ind_r_6 = find(sat_num==6);
ind_r_5 = find(sat_num==5);

error1 = Result(:,12:14);
if mean(error1(:,1)) < 0.01 && mean(error1(:,2)) < 0.01   
    en1 = zeros( length(error1(:,1)),1 );
    bias1 = 0;%mean(error1(:,1));
    bias2 = 0;%mean(error1(:,2));
    for i = 1:length(error1(:,1))
        en1(i) = norm([error1(i,1)-bias1, error1(i,2)-bias2]);
    end
    figure(12);
%     subplot(1,2,1)    
%     h = histfit(error1(:,1)-bias1);
    %delete(h(1));
    [nelements, centers] = hist(error1(:,1)-bias1, 30);
    grid on;
%     subplot(1,2,2)
%     h = histfit(error1(:,2)-mean(error1(:,2)));
%     grid on;
%     subplot(1,2,2)
    h = histfit(en1);
    grid on;
    delete(h(2));
    
    figure(14)
    [nelements, centers] = hist(en1, 30);
    bar(centers*1000, nelements./sum(nelements)*100, 1);
    grid on;
    xlabel('Average Horizontal Error (mm)');
    ylabel(['Percentage of total ', num2str(sum(nelements)),' trials.'])
else
    disp('Positioning error is too large!!!')
end