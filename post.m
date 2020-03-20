%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%    matlab code for post processing of HPC coursework
%%%%%%%%%%55

clc
clear
close all

Nx = 10;
Ny = 30;
% filename = ;
% opts.Sheet = '2007';
% opts.SelectedVariableNames = [1:3]; 
% opts.DataRange = '2:5';
% M = readmatrix('example.txt',opts);
M = load('streamFunctiom.txt');
N = load('vorticity.txt');
BC = load('vorticity_bc.txt'); %%left-right-bottom-top

Xs = reshape(M(:,1),[Nx-2,Ny-2]);
Ys = reshape(M(:,2),[Nx-2,Ny-2]);
S = reshape(M(:,3),[Nx-2,Ny-2]);
%Xv = reshape(N(:,1),[Nx-2,Ny-2]);
%Yv = reshape(N(:,2),[Nx-2,Ny-2]);
V = reshape(N(:,3),[Nx-2,Ny-2]);

X = zeros(Nx,Ny);
X(2:Nx-1,2:Ny-1) = Xs(:,:);
X(2:Nx-1,1) = Xs(:,1);
X(2:Nx-1,Ny) = Xs(:,Ny-2);
X(Nx,:) = X(Nx-1,:)+X(2,:);

Y = zeros(Nx,Ny);
Y(2:Nx-1,2:Ny-1) = Ys(:,:);
Y(1,2:Ny-1) = Ys(1,:);
Y(Nx,2:Ny-1) = Ys(1,:);
Y(:,Ny) = Y(:,2)+Y(:,Ny-1);

S_plot = zeros(Nx,Ny);
S_plot(2:Nx-1,2:Ny-1) = S(:,:);

 V_plot = zeros(Nx,Ny);
 V_plot(2:Nx-1,2:Ny-1) = V(:,:);
 V_plot(1,2:Ny-1) = BC(1:Ny-2);     %left
 V_plot(Nx,2:Ny-1) = BC((Ny-2)+1:2*(Ny-2));        %right
 V_plot(2:Nx-1,1) = BC(2*(Ny-2)+1:2*(Ny-2)+(Nx-2));            %bottom
 V_plot(2:Nx-1,Ny) = BC(2*(Ny-2)+(Nx-2)+1:length(BC));            %top

figure(1)
contourf(X,Y,S_plot);
% contourf(Xs,Ys,S,50, 'edgecolor','none');
% colormap jet
colorbar

figure(2)
contourf(Xs,Ys,V);
%contourf(X,Y,V_plot,50, 'edgecolor','none');colormap jet
colorbar;



%%for debug

%%% initialization matrix A
n = 5;
dx = 1/6;
dy = 1/6;
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
        
        
        s_in = [1.94791e-11    2.23274e-11    4.18686e-11    -7.91427e-11    -8.63992e-05    1.66306e-11    4.29768e-21    -5.18892e-21    1.5552e-07    -8.63969e-05    1.90805e-11    5.12641e-21    -1.69407e-21    1.5552e-07    -8.63995e-05    1.66303e-11    6.01385e-21    -3.45217e-21    1.5552e-07    -8.64021e-05    1.94795e-11    2.23287e-11    4.1873e-11    -7.91627e-11    -8.63992e-05   ] ;
        term1 = A\s_in';
        term1 = reshape(term1 , [n,n]);