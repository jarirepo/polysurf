% ************************************************************************
% polysurfdemo_public.m -- polysurf demo script
%
% Performs a bilinear interpolation between 4 connected 3-d parametric 
% functions. The functions are first discretized by sampling into polylines 
% and then interpolated using the Coon's bilinear interpolation method 
% which is implemented for polyline objects in the accompanying script 
% "polysurf.m". Further instructions can be found throughout the script
% which contains 3 primary steps. Have fun!
% 
% Dependencies:
%   polysurf.m
%   cpolysurf.mexw64
% 
% Date   : 10/2016
% Author : Jari Repo, University West, jari.repo@hv.se
% ************************************************************************
clear all, close all, clc

% ************************************************************************
% Step 1. Define the four boundary curves to the bilinear surface G(u,v) in
%         the XYZ-space. These are called G(u,0),G(u,1),G(0,v),G(1,v),
%         where u,v are the surface parameters in the uv-space with values
%         0<=u,v<=1
% 
% The following rules apply when defining the boundary edge functions:
% * curve G(u,0): goes in the positive u-direction at v=0, starts at G(0,0)
% and ends at point G(1,0).
% * curve G(u,1): goes in the positive u-direction at v=1, starts at G(0,1)
% and ends at point G(1,1).
% * curve G(0,v): goes in the positive v-direction at u=0, starts at G(0,0)
% and ends at point G(0,1).
% * curve G(1,v): goes in the positive v-direction at u=1, starts at G(1,0)
% and ends at point G(1,1).
% 
% Verify in the figure named "Definition curves" that the curves are 
% properly defined

% General settings:
% 
% use_mex           - set to 1 to call the MEX-implementation of polysurf,
%                     requires "cpolysurf.mexw64"
% plot_contour      - set to 1 to show the surface contour plot in Step 3
% surf_transp       - set transparency of the plotted surface using surf()
% surf_linestyle    - linestyle of output surface, try '-' to view the underlying meshing
% surfu             - number of segm. in the u-direction
% surfv             - number of segm. in the v-direction
% tighten           - set to include polyline points on the surface boundaries

use_mex = 0;
plot_contour = 1;
surf_transp = .95;
surf_linestyle = 'none';
surfu = 40;
surfv = 40;
tighten = 1;

% Here are 2 examples for you to explore
if 1
    % Example #1    
    % G(u,0) -- parabolically twisted abs. sine
    Gu0 = { @(t) t; ...
            @(t) -4*.5*t.*(t-1); ...
            @(t) abs(.25*sin(4*pi*t)) };

    % G(u,1) -- squared quarter sine
    Gu1 = { @(t) t; ...
            @(t) ones(size(t)); ...
            @(t) sin(pi/2*t).^2 };

    % G(0,v) -- circular arc
    G0v = { @(t) zeros(size(t)); ...
            @(t) (1/2)*(1-cos(pi*t)); ...
            @(t) sin(pi*t) };

    % G(1,v) -- parabolic
    G1v = { @(t) ones(size(t)); ...
            @(t) t; ...
            @(t) t.^2 };
else
    % Example #2
    % G(u,0) -- parabolically twisted abs. sine
    Gu0 = { @(t) t; ...
            @(t) -4*.5*t.*(t-1); ...
            @(t) abs(.5*sin(1*pi*t)) };

    % G(u,1) -- parabolically twisted abs. sine
    Gu1 = { @(t) t; ...
            @(t) 1+4*.35*t.*(t-1); ...
            @(t) abs(.5*sin(2*pi*t)) };

    % G(0,v) -- circular arc
    G0v = { @(t) zeros(size(t)); ...
            @(t) (1/2)*(1-cos(pi*t)); ...
            @(t) sin(pi*t) };

    % G(1,v) -- shifted parabolic
    G1v = { @(t) ones(size(t)); ...
            @(t) t; ...
            @(t) -4*.5*t.*(t-1) };
end

% ************************************************************************
% Goto Step 2 to modify the number of sampling points for each boundary
% function G(u,0),G(u,1),G(0,v),G(1,v)
% ************************************************************************

N = 101;
t = linspace(0,1,N);

% View the boundary curves
figure('name','Definition curves')
subplot(221),set(gca,'fontsize',9)
    plot(t,Gu0{1}(t), t,Gu0{2}(t), t,Gu0{3}(t))
    box off
    xlabel('t'),title('G(u,0)')
    legend('x(t)','y(t)','z(t)','location','best')
subplot(222),set(gca,'fontsize',9)
    plot(t,Gu1{1}(t), t,Gu1{2}(t), t,Gu1{3}(t))
    box off
    xlabel('t'),title('G(u,1)')
    legend('x(t)','y(t)','z(t)','location','best')
subplot(223),set(gca,'fontsize',9)
    plot(t,G0v{1}(t), t,G0v{2}(t), t,G0v{3}(t))
    box off
    xlabel('t'),title('G(0,v)')
    legend('x(t)','y(t)','z(t)','location','best')
subplot(224),set(gca,'fontsize',9)
    plot(t,G1v{1}(t), t,G1v{2}(t), t,G1v{3}(t))
    box off
    xlabel('t'),title('G(1,v)')
    legend('x(t)','y(t)','z(t)','location','best')

% 3-d plot    
fh3d = figure('name','3-d plot');
    set(gcf,'renderer','opengl','doublebuffer','on')
    set(gca,'fontsize',9)
    X = [Gu0{1}(t); Gu1{1}(t); G0v{1}(t); G1v{1}(t)]';
    Y = [Gu0{2}(t); Gu1{2}(t); G0v{2}(t); G1v{2}(t)]';
    Z = [Gu0{3}(t); Gu1{3}(t); G0v{3}(t); G1v{3}(t)]';
    ph = plot3(X,Y,Z);
    view(3)
    grid on
    xlabel('x'),ylabel('y'),zlabel('z')
    title('Surface boundary curves')
    legend('G(u,0)','G(u,1)','G(0,v)','G(1,v)','location','eastoutside')
    cameratoolbar
clear t

pause

% ************************************************************************
% Step 2. Sample the definition curves to convert them into polylines
%         (or sequences of points). Store them into a cell array
% 
%         Different number of sampling points can be used for the indivudal
%         surface boundary curves as shown in this example. 
%         Sufficient number of points is needed to represent the underlying
%         boundary definition functions.

t1 = linspace(0,1, 40);
t2 = linspace(0,1, 25);
t3 = linspace(0,1, 20);
t4 = linspace(0,1, 20);

% Sample and store the boundary functions into a 4-cell array
% Do not modify this!
P = cell(1,4);
P{1} = [Gu0{1}(t1); Gu0{2}(t1); Gu0{3}(t1)];    % G(u,0)
P{2} = [Gu1{1}(t2); Gu1{2}(t2); Gu1{3}(t2)];    % G(u,1)
P{3} = [G0v{1}(t3); G0v{2}(t3); G0v{3}(t3)];    % G(0,v)
P{4} = [G1v{1}(t4); G1v{2}(t4); G1v{3}(t4)];    % G(1,v)

% To view the content in P
% P{1},P{2},P{3},P{4}

% Output the individual polylines in the 3-d view
set(ph,'Visible','off') % hide the analytical functions
hold on
for k=1:4
    plot3(P{k}(1,:),P{k}(2,:),P{k}(3,:),'.:','Color',get(ph(k),'Color'))
end
hold off
title('Sampled surface boundary curves')

pause

% ************************************************************************
% Step 3. Generate a polysurf which is a bilinear interpolation between all
%         four polylines (piecewise linear functions)

disp('Generating surface ...')
if use_mex
    % try out the MEX implementation "cpolysurf.mexw64"
    % for 64-bit Windows systems
    [Gx,Gy,Gz] = cpolysurf(P,surfu,surfv,tighten); 
else
    [Gx,Gy,Gz] = polysurf(P,surfu,surfv,tighten);
end

[nv nu] = size(Gx);

hold on
% plot3(Gx,Gy,Gz,'.')
sh = surf(Gx,Gy,Gz);
set(sh,'linestyle',surf_linestyle,'facealpha',surf_transp)
if plot_contour
    % show contour plot on the "ground"
    [c,h] = contour(Gx,Gy,Gz,'LevelList',(0:.1:1));
end
hold off
axis equal
title(sprintf('Bilinear surface (%d faces)',(nu-1)*(nv-1)))

% add some shading
lighting phong
light

% --- End of demo program ---
