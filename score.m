function [score_tb] = score(mu,P,obsv,distance_tb)
% Lyuzhou Zhuang
% This function is to transform a distance table into a score table for data
% association functions
% score = distance*(1-probablity(distance)) probability~N(X,P)
% P = covariance matrix
    %sigma = sqrt(P);
    [m,n] = size(distance_tb);
    score_tb = NaN(m,n);
    density_tb = NaN(m,n); % a temporary table to store values for debugging
    for i = 1:m
        for j = 1:n
            density = mvnpdf(obsv(:,i)',mu(1:3,j)',P(1:3,1:3,j));
            density_tb(i,j) = density; % a temporary table to store values for debugging
            score_tb(i,j) = distance_tb(i,j)*(1-density);
        end
    end
end % add a breakpoint here for debugging