function [ obsv_labeled,scores ] = greedy( obsv,statemeans,cov )
% This code is greedy algorithm for associating observed points to the statemeans  
% Input arguments:-p = observed points matrix m=statemeans array; mean_ind = associated points index,
% points_i = observed points index, cov-statecovariances
% Output arguments:- ap=associated observed points, score = score of the matrix.

    [~,m] = size(obsv);                    %% finding size of observed points 
    [~,n] = size(statemeans);              %% finding size of statemeans
    distance_tb = zeros(m,n);             %% initializing the matrix
 
    %% Generate matrix for using in greedy algorithm
    % outer loop is the observed points loop and inner loop is the statemeans
    % loop. Therefore, we calculate distance of each statemeans from the
    % observed points.
    for i = 1:m      
        for j = 1:n
            tmp = obsv(:,i)-statemeans(1:3,j);    %% calculating the difference between mean and points
            distance_tb(i,j) = sqrt(tmp'*tmp);       %% calculating euclidean distance
        end
    end
    [R] = score(statemeans,cov,obsv,distance_tb);   %% calculating the score table which is used for the following algorithm
    %% Performing greedy algorithm
 
    R1 = R;             %% creating temporary matrix for performing greedy algorithm
    ind = [];           %% initializing index array which will contain index of min valued element in each rows
    scores = 0;     	%% initializing score
    for i = 1:m
        [a,b] = min(R1(i,:)); 	%% finding minimum value in each row(a=min value,b=position of the value in row)
        ind = [ind b];        	%% storing the position in
        R1(i,:) = NaN;        	%% deleting the row
        R1(:,b) = NaN;        	%% deleting the column with minimum value from row
        scores = scores+a;    	%% calculating the score (adding all the minimum values extracted from rows)
    end

    %% sorting the minimum index array in increasing order alongwith sorting of observed points index

    % here the minimum array index will already be associated to the observed
    % points as we are searching for the minimum values from each row(i.e distance of statemeans from each points)

    points_i = 1:m ;     	%% cretaing array of observed points index which is used for association
    mean_ind = ind;          	%% storing the indices of minimum values in a new array
    for j = 1:length(points_i)
        for i = 1:length(points_i)-1
            % sorting check is based on minimum indexes found
            if mean_ind(i) >= mean_ind(i+1)
                tmpp = mean_ind(i);
                mean_ind(i) = mean_ind(i+1);
                mean_ind(i+1) = tmpp;
                % sort the observed points index as well 
                tmpp = points_i(i);
                points_i(i) = points_i(i+1);
                points_i(i+1) = tmpp;
            end
        end    
    end

    %% association of observed points to statemeans(labels)
    % ap - Associated Points matrix. It will consist of associated observed points for every label 
    obsv_labeled=statemeans(1:3,:);                              %% initialize the matrix by positions of statemeans(1st 3 rows,i.e x,y,z coordinate)
    for i=1:length(mean_ind)
        obsv_labeled(:,mean_ind(i)) = obsv(:,points_i(i));   %% for every minimum index present for points substitute observed points values in matrix 
    end

end

