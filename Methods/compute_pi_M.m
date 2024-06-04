%Created on 2022.
%author:Louise Xiang

clear all
close all

path = '/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/catalog_all/';

%% Parameters from previous synthetic tests:
beta=log(10);    % valeur de b de la G-R * log(10)
m0=-5; % plus petite magnitude possible
m1=5; m2=6; % magnitudes min & max des "gros" chocs principaux (avec K variable). Les séismes m<m1 sont "normaux"%

%% 1- Load data:
load ('/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/catalog_all/amatrice.mat');
t = t - floor(t(1));

%% 2 - Cut-off at m0: 
J = find(ml >= m0);
t = t(J); ml = ml(J); x = lat(J); y = lon(J); z = z(J);

clear J

%% 3 - Compute smoothed magnitude mliss:
L = length(t);
dn = 100;
M = movmean(ml,dn);
I = find(ml > m1);

%% 4 - pi: Earthquake detection probability

%Fixed parameter:
euler_cst = 0.5772;

%step1 : Get smoothed magnitude (<=> Mmodele) 
mbar = M';

%Saving mbar used to compute pi later: 
baseFileName = 'M.dat';
fullFileName = fullfile(path, baseFileName);
fid=fopen(fullFileName,'w');
for n=1:length(mbar)
    fprintf(fid,'%d\n',mbar(n));
end
fclose(fid);

%step2: Compute mu, delta, pi:
mu = mbar - (euler_cst/beta);
delta0 = exp(-beta*(m0 - mu')); %%%%%%%% NOUVEAU
pi = (1./delta0).*(1- exp(-delta0)); %%%%%%%% NOUVEAU

%Save EQ dataframe with their associated pi:
baseFileName = 'data.dat';
fullFileName = fullfile(path, baseFileName);
fid=fopen(fullFileName,'w');
for n=1:L
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%e\n',t(n),ml(n),x(n),y(n),z(n),pi(n));   
end
fclose(fid);
clear n baseFileName fullFileName

%Save EQ with ml > m1:
baseFileName = 'list.dat';
fullFileName = fullfile(path, baseFileName);
fid=fopen(fullFileName,'w');
for n=1:length(I)
    fprintf(fid,'%d\n',I(n));
end
fclose(fid);

clear n baseFileName fullFileName
