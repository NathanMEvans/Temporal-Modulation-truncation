% Open ended dynamic time warping based on code by:
    % Paolo Tormene, Toni Giorgino, Silvana Quaglini, Mario Stefanelli
% with the license
    % Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
    % Signal Analysis and Machine Perception Laboratory,
    % Department of Electrical, Computer, and Systems Engineering,
    % Rensselaer Polytechnic Institute, Troy, NY 12180, USA

% Taking logic from 
% Matching incomplete time series with dynamic time warping: an algorithm and an application to post-stroke rehabilitation
% Paolo Tormene et al

% Dynamic time warping of two signals
% Modified cost funcion to inactive channels

% Nathan Evans, Newcastle University

function [distance_openEnded,fraction,truncation_point,distance]=dtw_openEnded(s,t,w)
% s: signal 1, size is ns*k, row for time, colume for channel 
% t: signal 2, size is nt*k, row for time, colume for channel 
% w: window parameter
%      if s(i) is matched with t(j) then |i-j|<=w
% distance_openEnded: resulting distance using open ended algorithm
% fraction:fraction used in open ended algorithm
% distance: resulting distance entire signal (regular DTW algorithm)

if nargin<3
    w=Inf;
end

ns=size(s,1);
nt=size(t,1);
if size(s,2)~=size(t,2)
    error('Error in dtw(): the dimensions of the two input signals do not match.');
end
w=max(w, abs(ns-nt)); % adapt window size

%% initialization
D=zeros(ns+1,nt+1)+Inf; % cache matrix
D(1,1)=0;

%% begin dynamic programming
for i=1:ns
    for j=max(i-w,1):min(i+w,nt)
        difference = (s(i,:)-t(j,:));
        active = any([s(i,:);t(j,:)]);
        % not changing the code now but I don't think difference(active) is
        % any different to difference
        % I also think for binary data maybe doing
        % sum(abs(difference))/sum(active) might make more sense
        
        %% USED FOR GEORGE'S ANALYSIS
        %cost = norm(difference(active))/sqrt(sum(active));

        %% FIXED /0 error
        if sum(active) == 0
            cost = 0;
        else
            cost = norm(difference(active))/sqrt(sum(active));
        end

        %% Ideal, don't know if I wanna use it
        % if sum(active) == 0
        %     cost = 0;
        % else
        %     cost = sum(abs(difference(active)))/sum(active);
        % end
        %cost=norm(s(i,:)-t(j,:));
        D(i+1,j+1)=cost+min( [D(i,j+1), D(i+1,j), D(i,j)] );
        
    end
end

v = (1:ns)./ns;
[distance_openEnded, j] = min(D(2:end, nt+1)./v');
fraction = (100*j/ns);
truncation_point = j;
distance=D(ns+1,nt+1);
