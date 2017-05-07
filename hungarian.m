function [ obsv_labeled,scores ] = hungarian( obsv,stmns,cov )
%This code does association of observed points to labels(statemeans) based
%on the hungarian method. Input arguments:- p - observed
%points,m-statemeans,cov-statecovariances Output arguments:- ap- associated points(labeled)
%matrix, score- score of the hungarian method

    %% step 1 creating R matrix
    [~,m] = size(obsv);                    %% finding size of observed points 
    [~,n] = size(stmns);                    %% finding size of statemeans
      
    distance_tb=zeros(m,n);                    %% initializing the matrix
    % outer loop is the observed points loop and inner loop is the statemeans
    % loop. Therefore, we calculate distance of each statemeans from the
 	% observed points.
    for i = 1:m      
        for j = 1:n
            co=obsv(:,i)-stmns(1:3,j);      %% calculating the difference between mean and points
            distance_tb(i,j)=sqrt(co'*co);       %% calculating euclidean distance
        end
    end
    [R] = score(stmns,cov,obsv,distance_tb);   %% calculating the score table which is used for the following algorithm
    %% Ensure that the matrix is square by the addition of dummy rows/columns if necessary. 
    [r0,c0]=size(R);                    %% finding size of R matrix
    R0=R;                               %% storing the matrix R
    % making R a square matrix by substituting rows and columns of identity
    % into unbalances matrix. Unbalanced matrix is created due to missing
    % observations.
    if r0<c0
        new=c0-r0;
        inew=eye(new,c0);
        R=[R; inew];
    end
    if c0<r0 
        new=r0-c0;
        inew=eye(r0,new);
        R=[R inew];
    end

    %% step 3 Reduce the rows by subtracting the minimum value of each row from that row.
    R1=R;                           %% using temporary matrix for editing
    [r,c]=size(R1);
    for i =1:r
        a=min(R1(i,:));                 %% finding minimum value in every row
        R1(i,:)=R1(i,:)-a;              %% subtracting the value from that particular row
    end

    %% step 4 If there are columns without a zero, reduce the columns by subtracting the minimum value of each column from that column
    R2=R1;                          %% using temporary matrix for column reduction
    for j =1:c
        a=min(R2(:,j));             %% finding minimum value in each column
        R2(:,j)=R2(:,j)-a;          %% subtracting it from that particular column
    end
    R3=R2;                          %% creating reduced matrix

    %% step 5 Cover the zero elements with the minimum number of lines it is possible to cover them with.
    asc=0;                          %% initailizing associations variable
    % if associations is not equal to the number of rows of matrix then keep
    % performing hungarian method
    while asc ~= r
        R4=R3;asc=0;                    %% initializing associations for the inner loop
        rscan=zeros(r,1);               %% row scan matrix which consisits of number of zeros in every row
        poszinrow=zeros(r,1);           %% this array contains poition of zeros in rows with only 1 zero
        for i=1:r
            if c-(length(find(R4(i,:)))) == 1 %% this loop is for the rows with only 1 zero
                rscan(i)=i;
                for j=1:c
                    if R4(i,j) == 0
                        poszinrow(i) = j;       %% position of zeros by scanning columns of that row
                        asc=asc+1;              %% if zero found, associate it with the row and update associations
                        for k=1:r
                            if k~=i 
                                R4(k,j) = NaN;  %% delete the other column elements, except the zero
                            end
                        end
                    end
                end
            end
        end
        R5=R4;                                  %% save the result in temporary matrix
        iR5=sum(isnan(R5));                     %% matrix that contains position of deleted elements
        cscan=zeros(1,c);                       %% col scan matrix which consisits of number of zeros in every col
        poszincol=zeros(1,c);                   %% this array contains poition of zeros in col with only 1 zero
        for j=1:c
            if iR5(j) == 0                      %% this loop is for the columns with only 1 zero
                cscan(j)=j;
                for i=1:r
                    if R5(i,j) == 0
                        poszincol(j) = i;       %% position of zeros by scanning rows of that column
                        asc=asc+1;              % if zero found, associate it with the col and update associations
                        for l=1:c
                            if l~=j 
                                R5(i,l) = NaN;  %% delete the other row elements, except the zero
                            end
                        end
                    end
                end
            end
        end

        %% step 6 Add the minimum uncovered element to every covered element.
        % Perform above row and column deletion to generate matrix which will be
        % used for comparison, but this matrix will contain all the elements that
        % will deleted for the selected associations rather than leaving just a
        % zero undeleted like the previous matrix
        R5r=R3;                                 %% save the original result after reductions into another temporary matrix
        for i=1:r
            if c-(length(find(R5r(i,:)))) == 1  %% this loop is for the rows with only 1 zero
                for j=1:c
                    if R5r(i,j) == 0
                        for k=1:r
                            R5r(k,j) = NaN;    %% delete the other column elements alongwith the zero
                        end
                    end
                end
            end
        end      
        iR5r=sum(isnan(R5r));                   %% matrix that contains position of deleted elements
        for j=1:c
            if iR5r(j) == 0                     %% this loop is for the columns with only 1 zero
                for i=1:r
                    if R5r(i,j) == 0
                        for l=1:c
                            R5r(i,l) = NaN;     %% delete the other row elements alongwith the zero
                        end
                    end
                end
            end
        end
        mvR5=min(min(R5r));                     %% the above matrix gives the uncovered elements from which we find the minimum value
        R6=R3;                                  %% Use R6 matrix for performing editing operations for step 6 of hungarian method
        Rn=R5r-mvR5;                            %% Subtract minimum uncovered value from temporary matrix and store it in reference matrix
        iRn=isnan(Rn);                          %% find the location of deleted elements from reference matrix
        % adding minimum value of open elements to intersecting points. This is
        % done by finding the location of vertical and horizontal lines.
        vv = find( poszinrow );     % finding location of vertical lines(deleted columns)
        hv = find( poszincol );     % finding location of horizontal lines(deleted rows)
        % Adding minimum  value to the intersecting points
        for i =1:length(hv)
            for j=1:length(vv)
                R6(poszincol(hv(i)),poszinrow(vv(j)))=R6(poszincol(hv(i)),poszinrow(vv(j))) + mvR5; 
            end
        end
        % generating new matrix for performing steps from step 5 again if number of
        % associations is not equal to the number of rows.
        % Substituting rest of the elements from the covered(deleted columns 
        % and rows) as it is in the reference matrix
        for i=1:r
            for j=1:c
                if iRn(i,j)==1
                    Rn(i,j)=R6(i,j);        
                end
            end
        end
        R3=Rn;                          % Make reference matrix as the new original matrix that will start from step 5
    end
    
    % If the associations are complete start assigning the observed points to
    % labels(statemeans)
    asm=zeros(size(R1));        %% initailize assignment matrix,by making cell values in matrix 0(off)
    for i=1:length(poszinrow)
        if poszinrow(i) ~= 0    %% use poszinrow(position in row) matrix to find location of zero
            asm(i,poszinrow(i))=1;  %% for the particular row if zero is assigned then make the cell 1(on) in the particular column
        end
    end

    for j=1:length(poszincol)
        if poszincol(j) ~=0     %% use poszinrow(position in col) matrix to find location of zero
            asm(poszincol(j),j)=1;  %% for the particular column if zero is assigned then make the cell 1(on) in the particular row
        end
    end

    % Sometimes the associations will be complete but associated element will
    % not be calculated. Hence,check for empty assignment. 
    uar=[];                     %% initialize unassigned row matrix
    for i=1:r
        if sum(asm(i,:)) == 1   %% this is scanning through row to find empty assignment
            continue
        else
            uar=[uar i];        %% if found store the row number
        end
    end
    uac=[];                     %% initialize unassigned col matrix
    for j=1:c
        if sum(asm(:,j)) == 1   %% this is scanning through column to find empty assignmen
            continue
        else
            uac=[uac i];        %% if found store the column number
        end
    end

    if (isempty(uar) == 0) && (isempty(uac) == 0)
    iter=1;                     %% temporary variable for column assignment indexing
        for i=1:length(uar)
            asm(uar(i),uac(iter))=1;    %% in the assignment matrix, assign the unassigned column to unassigned row.
            iter=iter+1;            %% increase the index for next assignment
        end
    end

    %% assigning means to proper points and calculating score
    scores=0;                    %% this will be score, which is used for qualitative analysis
    [ra,ca]=size(asm);          %% find the size of assignment matrix
    obsv_labeled=stmns(1:3,:);              %% initialize the associated points matrix with 1st three rows of statemeans matrix
    % the below is an if statement because the original matrix(score table) which is used
    % for hungarian method is rectangular(unbalanced) due to missing
    % observation and it was made square.
    if r0<ra
        for j=1:ca
            for i=1:r0
                if asm(i,j)==1  %% if assigned zero is present in assignmnet matrix then assign the particular observed points to statemeans 
                    obsv_labeled(:,j)=obsv(:,i);     %% particular observed points is substitued for the particular statemeans
                    scores=scores+R0(i,j);    %% also the value from the score table is added for calculating score
                end
            end
        end
    end
    if c0<ca
        for j=1:c0
            for i=1:ra
                if asm(i,j)==1  %% if assigned zero is present in assignmnet matrix then assign the particular observed points to statemeans 
                    obsv_labeled(:,j)=obsv(:,i);     %% particular observed points is substitued for the particular statemeans
                    scores=scores+R0(i,j);    %% also the value from the score table is added for calculating score
                end
            end
        end
    end    
    if (r0 == ra) && (c0 == ca)
        for j=1:ca
            for i=1:ra
                if asm(i,j)==1  %% if assigned zero is present in assignmnet matrix then assign the particular observed points to statemeans 
                    obsv_labeled(:,j)=obsv(:,i);     %% particular observed points is substitued for the particular statemeans
                    scores=scores+R0(i,j);        %% also the value from the score table is added for calculating score
                end
            end
        end
    end

end

