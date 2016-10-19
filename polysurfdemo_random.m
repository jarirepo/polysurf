% polysurfdemo_random.m
% 
% Dependencies: polysurf.m
% 
% Date   : 10/2016
% Author : Jari Repo, University West, jari.repo@hv.se
% ************************************************************************
clear all, close all, clc

% Specify surface resolution in the u- and v-directions
surfu = 40;
surfv = surfu;
transp = 1;
N = 19;
A = .5;         % max. Z-level
plot_contour = 1;

% Generate randomized point sequences representing the surface boundary
% curves
% Gu0 = zeros(3,N);
% Gu1 = zeros(3,N);
% G0v = zeros(3,N);
% G1v = zeros(3,N);

t = linspace(0,1,N);
f = 4*A*t.*(t-1).^2;
Gu0 = [t; zeros(1,N); (1/2)*(1+rand(1,N)).*f];   % G(u,0)
Gu1 = [t; ones(1,N); (1/2)*(1+rand(1,N)).*f];    % G(u,1)
G0v = [zeros(1,N); t; (1/2)*(1+rand(1,N)).*f];   % G(0,v)
G1v = [ones(1,N); t; (1/2)*(1+rand(1,N)).*f];    % G(1,v)

% Define surface corners by averaging
P00 = (Gu0(:,1)+G0v(:,1))/2;
P10 = (Gu0(:,end)+G1v(:,1))/2;
P01 = (G0v(:,end)+Gu1(:,1))/2;
P11 = (Gu1(:,end)+G1v(:,end))/2;
% Ensure endpoint connectivity
Gu0(:,1) = P00;
G0v(:,1) = P00;
Gu0(:,end) = P10;
G1v(:,1) = P10;
G0v(:,end) = P01;
Gu1(:,1) = P01;
Gu1(:,end) = P11;
G1v(:,end) = P11;

% Store the piecewise linear boundary functions into a 4-cell array
% This also sets the 3rd dimension to zero
P = cell(1,4);
P{1} = Gu0; % G(u,0)
P{2} = Gu1; % G(u,1)
P{3} = G0v; % G(0,v)
P{4} = G1v; % G(1,v)

% To view the content in P
% P{1},P{2},P{3},P{4}

% Generate the bilinear surface G(u,v)
[Gx,Gy,Gz] = polysurf(P,surfu,surfv,1);

% 3-d plot    
fh3d = figure('name','3-d plot');
    set(gcf,'renderer','opengl','doublebuffer','on')
    set(gca,'fontsize',9)
    % Output the individual polylines in the 3-d view
    hold on
    for k=1:4
        plot3(P{k}(1,:),P{k}(2,:),P{k}(3,:),'.:')
    end
    
    % plot3(Gx,Gy,Gz,'.')
    sh = surf(Gx,Gy,Gz);
    set(sh,'linestyle','none','facealpha',transp)
    if plot_contour
        % show contour plot on the "ground"
        [c,h] = contour(Gx,Gy,Gz,'LevelList',(0:.1:1));
    end
    hold off
    grid on
    view(3)
    axis equal
    cameratoolbar
    xlabel('x'),ylabel('y'),zlabel('z')
%     legend('G(u,0)','G(u,1)','G(0,v)','G(1,v)',...
%         'location','eastoutside')

    % add some shading
    lighting phong
    light
    zlim([0 1.1*A])
    title('Bilinear surface based on randomized boundaries')
    
    