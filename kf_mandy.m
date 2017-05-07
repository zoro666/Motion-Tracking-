function [ X_updated_prediction,P_updated_prediction ] = kf_mandy( X,P,obsv )    
    %% Initialization 
    % X = statemeans P = covariances

    F= [1 0 0 1 0 0;
        0 1 0 0 1 0;
        0 0 1 0 0 1;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1];
    
    H= [1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0];
    
    Q = diag([.1 .1 .1 .1 .1 .1]);
    R = diag([.1 .1 .1]);
    R_increased = diag([.4 .4 .4]);
    
    %% predict

    P_predicted = zeros(size(P));
    X_predicted = zeros(size(X));
    [~,~,num_of_labels] = size(P);
    
    % predicted state (mean)
    X_predicted = F*X;

    % predicted estimate covariance
    for i = 1 : num_of_labels
        P_predicted(:,:,i) = (F*P(:,:,i)*F')+Q;
    end

    %% Update
    % innovation or measurement residue
    [~,cxp]=size(X_predicted);
    for i = 1:cxp
        y(:,i)=obsv(:,i)-(H*X_predicted(:,i));
    end

    % for kalman gain
    Sk = zeros(3,3,num_of_labels);
    Kk = zeros(6,3,num_of_labels);
    if sum(y) == 0
        for i = 1:num_of_labels
            Sk(:,:,i) = (H*P_predicted(:,:,i)*H')+R_increased; %%innovation (or residual) covariance
            Kk(:,:,i) = P_predicted(:,:,i)*H'*(inv(Sk(:,:,i))); %% kalman gain
        end
    else
        for i = 1:num_of_labels
            Sk(:,:,i) = (H*P_predicted(:,:,i)*H')+R; %%innovation (or residual) covariance
            Kk(:,:,i) = P_predicted(:,:,i)*H'*(inv(Sk(:,:,i))); %% kalman gain
        end
    end

    % for updating X & P predition
    X_updated_prediction = zeros(size(X_predicted));
    P_updated_prediction = zeros(size(P));
    I = eye(6,6);
    for i = 1:num_of_labels
        X_updated_prediction(:,i) = X_predicted(:,i) + Kk(:,:,i)*y(:,i);
        P_updated_prediction(:,:,i) = (I - (Kk(:,:,i)*H))*P_predicted(:,:,i);
    end
end

