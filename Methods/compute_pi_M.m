clear all
close all

path = '/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/catalog_all/';

%% paramètres issus de "synthetique_modif_v4.m" 
beta=log(10);    % valeur de b de la G-R * log(10)
% p=1.1; c=10^(-5); % paramètres loi d'omori
% alpha=2; % paramètres de productivité
% gamma=2; % paramètre du kernel en distance (exposant de décroissance avec r)
% L0=0.5*10^(-2);     % rayon de rupture (en km) pour un magnitude = 0 
% tmax=365;
% Lmax=20; Zmax=10; % les + gros séismes (magnitudes entre m1 et m2) sont dans une boite [0 Lmax] x [0 Lmax] x [0 Zmax] (en km)
m0=-5; % plus petite magnitude possible
% k=0.001*exp(-alpha*m0);
m1=5; m2=6; % magnitudes min & max des "gros" chocs principaux (avec K variable). Les séismes m<m1 sont "normaux"%
% b=beta/(alpha-beta)*k*exp(beta*m0)/(1-exp(-beta*(m1-m0)))*(exp((alpha-beta)*m1)-exp((alpha-beta)*m0))*c^(1-p)/(p-1);
% s=sprintf('taux de branchement estime = %f (probleme si >1 -> diminuer la valeur de k). ',b); disp(s);         % attention : si >1 alors divergence !

%% charger les données :
%
load ('/home/xianglo/Bureau/DATA_AMATRICE/DBSCAN_David/new_all_run_files_for_etas/catalog_all/amatrice.mat');

t = t - floor(t(1));

%% couper le catalogue à m0 : 
J = find(ml >= m0);

t = t(J); ml = ml(J); x = lat(J); y = lon(J); z = z(J);

clear J

%% 5 - mliss : magnitude lissée
L = length(t);

dn = 100;

% Utilisation de "movmean" :
M = movmean(ml,dn);

I = find(ml > m1);

%% 6 - pi : probabilité de détection

%paramètres choisis à partir de synthetique.m :
euler_cst = 0.5772;

%step1 : calcul de mbar lissée (<=> Mmodele) 
mbar = M';

%sauvegarder mbar utilisé pour calculer les nvelles valeurs de pi : 
baseFileName = 'M.dat';
fullFileName = fullfile(path, baseFileName);
fid=fopen(fullFileName,'w');
for n=1:length(mbar)
    fprintf(fid,'%d\n',mbar(n));
end
fclose(fid);

%step2 : calcul mu, delta, pi :
mu = mbar - (euler_cst/beta);
delta0 = exp(-beta*(m0 - mu')); %%%%%%%% NOUVEAU
pi = (1./delta0).*(1- exp(-delta0)); %%%%%%%% NOUVEAU

%sauvegarde du catalogue après calcul de proba de détection :
baseFileName = 'data.dat';
fullFileName = fullfile(path, baseFileName);
fid=fopen(fullFileName,'w');
for n=1:L
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%e\n',t(n),ml(n),x(n),y(n),z(n),pi(n));   
end
fclose(fid);
clear n baseFileName fullFileName

%Sauvegarde de la liste des ml > m1 séismes :
baseFileName = 'list.dat';
fullFileName = fullfile(path, baseFileName);
fid=fopen(fullFileName,'w');
for n=1:length(I)
    fprintf(fid,'%d\n',I(n));
end
fclose(fid);

clear n baseFileName fullFileName
