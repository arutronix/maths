
% MATLAB code to plot the MPI cannon matrix product results
% author: Aravind I.
% date: 11/13/2018
clear; clear all;
N = 4096;

%load the exec time into matlab workspace
%exec_time_no_transpose = [544.341190 55.494490 23.557356];
exec_time = [138.395898 34.640417 18.331599];
procs = [4 16 64];
% Calculation of FLOPS for matrix multiply: 2N / execution time
FLOPS = 2*N*N*N;
MFLOPS = (FLOPS./exec_time)/1e6;

% Calculate speed up
speedup = 558.210689./exec_time;

% Calculate parallel efficiency
parEff = 100*speedup./procs;

% Execution Time v. P
figure
plot(procs,exec_time,'-rs','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',5)
title('MPI Cannon Matrix Product: execution time v. procs')
xlabel('procs')
ylabel('execution time (s)')
grid on

% Plot Speed UP v. N
figure
plot(procs,speedup, '-rs','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',5)
title('MPI Cannon Matrix Product: speed up v. procs')
xlabel('procs')
ylabel('speed up')
grid on

% Plot Parallel Efficiency v. N
figure
plot(procs,parEff,'-rs','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',5)
title('MPI Cannon Matrix Product: parallel efficiency v. procs')
xlabel('procs')
ylabel('parallel efficiency (%)')
grid on

% Plot MFLOPS v. N
figure
plot(procs,MFLOPS,'-rs','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',5)
title('MPI Cannon Matrix Product: MFLOPS v. procs')
xlabel('procs')
ylabel('MFLOPS')
grid on

