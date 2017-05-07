% SampleCode.m  - for Project 1 in Advanced Computer Vision, Spring 17
%
% The goal of the project is to implement a Kalman Filter that tracks and
% labels 39 motion capture markers automatically.
%
% Given:
%	initc3d - a single frame that has correctly labeled markers
%	datac3d - sequence of frames of data that needs to be corrected
%
% Desired:
%	outputc3d - this is a data structure with the corrected marker labeling
%
% This sample code shows how to read and visualize the motion capture data
% and gives a rough outline of the overall tracking code framework.


% Clear workspace and initia  lize variables
close all; clear all; clc;
addpath(genpath('./')); %must add path w/subfolders for code to work
dir = '../08-21-15_Subject1/';
name = '24Form-Part1-Take1';

%% Initialization and Loading Data 
fps = 100; % # of frames to be processed
initframe = 827;
% Set data file paths
init_fn = [dir name '.frame_827.c3d']; % single frame w/ 39 correct labels
data_fn = [dir name '.labeled_auto.c3d']; % data to be corrected

% Load the input c3d files as matlab structures
initc3d = readMocapData39(init_fn); % load ground truth labeling of 1 frame 
initialvalues = initc3d.data;
datac3d = readMocapData39(data_fn); % load noisy observed data sequence
observations = datac3d.data;

% initialize output structure with observation sequence
outputc3d = datac3d;       
% however, wipe out all the marker locations (set them to NaN values)
outputc3d.data = nan(size(outputc3d.data));
% except for the ground truth frame marker values
outputc3d.data(initframe,:) = initialvalues;


% % Examples of how to generate visualizations
% % display the ground truth labeling frame
% p = visualizeSkeleton(initc3d,'init.avi');  
% % display data frame 827 
% p = visualizeSkeleton(datac3d,'data1frame.avi',initframe,initframe);  
% % display 3 seconds of data starting at frame 827
% fps = datac3d.freq;
% p = visualizeSkeleton(datac3d,'data3seconds.avi',initframe,initframe+3*fps);


% note: 3D marker points are stored in a big array, e.g. datac3d.data 
% each row represents the coords of the marker points in one frame.
% order of columns is x1,y1,z1,x2,y2,z2,x3,y3,z3,...,x39,y39,z39 so there
% are 3*39 = 117 columns in that data array.  Any points that are missing
% in the data have x,y,z values of NaN (which represents "Not a Number").
% numeric values are represented in millimeters.


% General sketch of tracking code

% Initialize mean values of state vectors from ground truth frame.
% columns are [x,y,z,dx,dy,dz] with location set from the labeled
% frame marker locations and the velocity set to zero.
statemeans = zeros(6,39); 
ptr = 1;
for i = 1:39
    statemeans(1:3,i) = initialvalues(1,ptr:(ptr+2))';
    ptr = ptr+3;
end

% Initialize covariance matrices for each state vector.  Each covariance
% is a 6x6 matrix.  The 39 of them are stored in a 6x6x39 array. 
statecovmats = zeros(6,6,39);
for i=1:39
    % the following is just an example.  You might find that something
    % else is more sensible.  This example sets initial covariance mat to
    % be diagonal, with location variances 10mm and velocity variances 1m.
    statecovmats(:,:,i) = diag([.9 .9 .9 .9 .9 .9]);
end

%%
% Loop through observations forward in time from initial frame
startframe = initframe+1;
endframe = startframe+fps;  %very short test, 1 second
score_tb = []; %% An array store scores of all frames processed

for currframe = startframe:endframe
    %% Load Observations
    % get valid observations from current frame into an array. the columns
    % will be [x,y,z]' observed locations. Note that number of observations
    % might be less than 39, because we are ignoring NaN values.
    observed = zeros(3,39);
    ptr = 1; colptr = 1;
    for i=1:39
        observed(1:3,colptr) = observations(currframe,ptr:(ptr+2))';
        ptr = ptr+3;
        if not(isnan(observed(1,colptr)))
            colptr = colptr+1;
        end
    end
    observed = observed(:,1:(colptr-1));  %shrink to size of valid data 
    
    %% Data Association (task_1)
    % Data Association code would go here, for determining which column
    % of observed data goes with which column of state vector array.  This
    % would be based on predicting each state forward one time step and
    % then associating predicted locations with observed locations.
    
    [observ_labeled,score] = greedy( observed,statemeans,statecovmats );
    score_tb = [score_tb score];  %% store scores of all frames processed
    
    %% Kalman Filter (task_2)
    % Kalman update code would go here, to update each state mean and
    % covariance with the associated observation determined above.
    % If there was no associated observation found, that state would
    % stay at its predicted location but with increased covariance to 
    % reflect additional uncertainty.  
    [ statemeans,statecovmats ] = kf_mandy( statemeans,statecovmats,observ_labeled );

    %% Save results
    % Finally, fill in the output marker values with the computed state
    % locations.  The maximum a posteriori location estimates are just
    % the location portions of the updated state mean vectors. Note that
    % we are killing several birds with one stone here: 1) labeling;
    % 2) smoothing; and 3) filling small gaps (using predicted locations), 
    % all by just setting output locations to computed mean locations.
    % save statemeans updated by the above code
    ptr = 1;
    for i=1:39
        %save statemeans updated by the above code
        outputc3d.data( currframe, ptr:(ptr+2) ) = statemeans(1:3,i)';
        ptr = ptr+3;
    end
end

%% Save and output results
% write out the "corrected" motion capture data
save 'corrected_marker_labels_greedy.mat' outputc3d

% visualize the corrected data
p = visualizeSkeleton(outputc3d,'datacorrected.avi',initframe,endframe);

% plot of the scores of motion capture data
plot(score_tb)
title('scores by greedy data association')
xlabel('frames')
ylabel('scores')
