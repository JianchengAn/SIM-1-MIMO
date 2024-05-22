function [ PA_WF ] = WF( P_total, sigma2, a )
%% WF: water-filling solution
% P_total: total transmit power
% sigma2: noise power
% a: singular value vector of MIMO channel;
% PA_WF: power allocation solution for each stream using WF;
P_sum = 0;
for ii = 1 : length(a) - 1
    P_sum = P_sum + ii*(sigma2/a(ii+1)^2 - sigma2/a(ii)^2); % Calculate the sum power when all (1:ii) channels are filled
    if P_sum >= P_total
        break
    end
end
if P_sum >= P_total
    PA_WF = sigma2/a(ii+1)^2 - sigma2./a(1:ii).^2 - (P_sum - P_total)/ii; % Subtract the excess power for all (1:ii) channels
    PA_WF = [PA_WF; zeros(length(a) - length(PA_WF), 1)]; % Match the dimension
else
    PA_WF = sigma2/a(ii+1)^2 - sigma2./a(1:ii+1).^2 + (P_total - P_sum)/(ii+1); % Add the extra power for all (1:ii+1) channels
end
end