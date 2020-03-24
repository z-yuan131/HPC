%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%    matlab code for post processing of HPC coursework
%%%%%%%%%%55

clc
clear
close all

Nx = 161;
Ny = 161;
% filename = ;
% opts.Sheet = '2007';
% opts.SelectedVariableNames = [1:3]; 
% opts.DataRange = '2:5';
% M = readmatrix('example.txt',opts);
M = load('streamFunctiom.txt');
N = load('vorticity.txt');
BC = load('vorticity_bc.txt'); %%left-right-bottom-top

Xs = reshape(M(:,1),[Ny-2,Nx-2]);
Ys = reshape(M(:,2),[Ny-2,Nx-2]);
S = reshape(M(:,3),[Ny-2,Nx-2]);
V = reshape(N(:,3),[Ny-2,Nx-2]);
Xs = Xs;
Ys = Ys;
S = S;
V = V;

X = zeros(Ny,Nx);
X(2:Ny-1,2:Nx-1) = Xs(:,:);
X(2:Ny-1,Nx) = Xs(:,Nx-2)+Xs(:,1);
X(1,:) = X(2,:);
X(Ny,:) = X(Ny-1,:);

Y = zeros(Ny,Nx);
Y(2:Ny-1,2:Nx-1) = Ys(:,:);
Y(Ny,:) = Y(2,:)+Y(Ny-1,:);
Y(:,Nx) = Y(:,Nx-1);
Y(:,1) = Y(:,Nx-1);

S_plot = zeros(Ny,Nx);
S_plot(2:Ny-1,2:Nx-1) = S(:,:);

 V_plot = zeros(Ny,Nx);
 V_plot(2:Ny-1,2:Nx-1) = V(:,:);
 V_plot(2:Ny-1,1) = BC(1:Ny-2);     %left
 V_plot(2:Ny-1,Nx) = BC((Ny-2)+1:2*(Ny-2));        %right
 V_plot(1,2:Nx-1) = BC(2*(Ny-2)+1:2*(Ny-2)+(Nx-2));            %bottom
 V_plot(Ny,2:Nx-1) = BC(2*(Ny-2)+(Nx-2)+1:length(BC));            %top

figure(1)
contourf(X,Y,S_plot);
% contourf(Xs,Ys,S,50, 'edgecolor','none');
% colormap jet
colorbar
title('stream');

figure(2)
contourf(Xs,Ys,V,50);
% contourf(X,Y,V_plot,50, 'edgecolor','none');colormap jet
colorbar;
title('vorticity')
