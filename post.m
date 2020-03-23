%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%    matlab code for post processing of HPC coursework
%%%%%%%%%%55

clc
clear
close all

Nx = 15;
Ny = 10;
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
%Xv = reshape(N(:,1),[Nx-2,Ny-2]);
%Yv = reshape(N(:,2),[Nx-2,Ny-2]);
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
% Y(1,2:Ny-1) = Ys(1,:);
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
contourf(Xs,Ys,V);
% contourf(X,Y,V_plot,50, 'edgecolor','none');colormap jet
colorbar;
title('vorticity')



%%for debug

%%% initialization matrix A
n = (Nx-2);
dx = 1/(Nx-1);
dy = 1/(Ny-1);
        B = zeros(n , n);
        A_ini = cell(n , n);
        A_ini(: , :) = {zeros(n)};

        for i = 1 : n
            for j = 1 : n
                if i == j
                    B(i , j) = 2.0/(dx*dx) + 2.0/(dy*dy); %A
%                     B(i , j) = 0; %B
%                     B(i , j) = 0; %C
                    if i < n
                        B(i + 1 , j) = -1.0/(dy*dy);
%                         B(i + 1 , j) = -1.0/2.0/dy;
%                         B(i + 1 , j) = 0;%defined in report 
                        B(i , j + 1) = -1.0/(dy*dy);
%                         B(i , j + 1) = 1.0/2.0/dy;
%                         B(i , j + 1) = 0;
                    end;
                end
            end;
        end;


        for i = 1 : n
            for j = 1 : n
                if i == j
                    A_ini{i , j} = B;
                    if i < n
                        A_ini{i + 1 , j} = -eye(n)*1.0/(dx*dx);     %defined in report    
%                         A_ini{i + 1 , j} = -eye(n)*0;     %defined in report  
%                         A_ini{i + 1 , j} = -eye(n)*1.0/2.0/dx;     %defined in report  
                        A_ini{i , j + 1} = -eye(n)*1.0/(dx*dx);
%                         A_ini{i , j + 1} = -eye(n)*0;
%                         A_ini{i , j + 1} = eye(n)*1.0/2.0/dx;
                    end;
                end
            end;
        end;

        A = cell2mat(A_ini);
        
        
        s_in = [-2.35145e-07    -5.00557e-07    -8.32109e-07    -1.27983e-06    -1.92527e-06    -2.92422e-06    -4.62662e-06    -7.95612e-06    -4.40051e-07    -9.35005e-07    -1.5481e-06    -2.36204e-06    -3.49714e-06    -5.14517e-06    -7.62592e-06    -1.10202e-05    -5.90084e-07    -1.25132e-06    -2.06324e-06    -3.12308e-06    -4.55598e-06    -6.52418e-06    -9.21546e-06    -1.23081e-05    -6.69003e-07    -1.41694e-06    -2.33045e-06    -3.51106e-06    -5.07942e-06    -7.17067e-06    -9.89822e-06    -1.28018e-05    -6.69027e-07    -1.417e-06    -2.33057e-06    -3.51129e-06    -5.07984e-06    -7.17141e-06    -9.89942e-06    -1.28034e-05    -5.90146e-07    -1.25147e-06    -2.06355e-06    -3.12368e-06    -4.55712e-06    -6.52626e-06    -9.21903e-06    -1.23135e-05    -4.40123e-07    -9.35183e-07    -1.54846e-06    -2.36277e-06    -3.49859e-06    -5.14801e-06    -7.63143e-06    -1.10302e-05    -2.35192e-07    -5.00675e-07    -8.32358e-07    -1.28035e-06    -1.92634e-06    -2.92654e-06    -4.63201e-06    -7.96982e-06  ] ;
        term1 = A*s_in';
        term1 = reshape(term1 , [n,n]);