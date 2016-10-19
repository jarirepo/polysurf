function [Gx,Gy,Gz] = polysurf(P, surfu, surfv, matchBoundaries)
%POLYSURF Bilinear interpolation between 4 polylines
% 
% USAGE: [Gx,Gy,Gz] = polysurf(P, surfu, surfv, [matchBoundaries=0])
% 
% INPUTS
% P         -- a 4-element cell array containing the definition points of
%              the four polylines. Each polyline is a sequence of points
%              that represent a boundary edge/curve G(u,v) of the final 
%              bilinear surface G(u,v) = [x(u,v),y(u,v),z(u,v)]
% 
%              The organization of P is as follows.
%              P{1}: curve G(u,0)
%              P{2}: curve G(u,1)
%              P{3}: curve G(0,v)
%              P{4}: curve G(1,v)
% 
% surfu     -- number of surface segments in the U-direction (16)
% surfv     -- number of surface segments in the V-direction (16)
% 
% matchBoundaries -- set to 1 to include the polyline control points in the
%                    generated surface. The number of segments in the u- 
%                    and v-directions may then become larger than specified 
%                    by surfu and surfv. On the other hand, the resolution
%                    will also be slightly improved near the definition
%                    points. Default is matchBoundaries=0
% 
% OUTPUTS
% Gx,Gy,Gz  -- matrices containing the X-,Y- and Z-values of the generated 
%              surface. These can be used as input to the built-in SURF()
%              function.

% Date   : 10/2016
% Author : Jari Repo, University West, jari.repo@hv.se
% 
% Modified 2016-10-19, JRE
% The boundary definition points can be included in the surface
% interpolation by setting matchBoundaries=1

if ~exist('surfu'),surfu = 16;end
if ~exist('surfv'),surfv = 16;end
if ~exist('matchBoundaries'),matchBoundaries = 0; end
if length(P)~=4,error('P must contain 4 polylines'),end

% Init outputs
Gx = []; Gy = []; Gz = [];

% Generated the u- and v-grids; 0<=(u,v)<=1
u = linspace(0,1,surfu+1);
v = linspace(0,1,surfv+1);

% Parametrize each polyline in the interval t=[0,1] based on the lengths of
% the individual segments
S(1:4) = struct('p',[],'t',[],'idx',[],'h',[]);
for k=1:4    
    S(k).p = P{k};
    S(k).t = parametrize( P{k} );    
end
    
if matchBoundaries
    % Inject the polyline parameters 0<t<1 into the u- and v-grid to
    % include the polyline definition points in the final surface.
    tol = 1e-6;
    for i=2:length(S(1).t)-1
        if ~any(abs(u-S(1).t(i)) < tol)
            u = [u, S(1).t(i)];
    end
    for i=2:length(S(2).t)-1
        if ~any(abs(u-S(2).t(i)) < tol)
            u = [u, S(2).t(i)];
        end
    end
    for i=2:length(S(3).t)-1
        if ~any(abs(v-S(3).t(i)) < tol)
            v = [v, S(3).t(i)];
        end
    end
    for i=2:length(S(4).t)-1
        if ~any(abs(v-S(4).t(i)) < tol)
            v = [v, S(4).t(i)];
        end
    end
    u = sort(u,'ascend');
    v = sort(v,'ascend');
end

% Map the parameter values (u,v) to boundary curve segment numbers to speed
% up later curve interpolation.
[S(1).idx, S(1).h] = mapParam(S(1).t, u);
[S(2).idx, S(2).h] = mapParam(S(2).t, u);
[S(3).idx, S(3).h] = mapParam(S(3).t, v);
[S(4).idx, S(4).h] = mapParam(S(4).t, v);

% Generate the bilinear surface
[Gx,Gy,Gz] = generateSurface(S,u,v);
end

% --- Helper functions ---
function [t] = parametrize(p)
    t = [0,cumsum(sqrt(sum((p(:,2:end)-p(:,1:end-1)).^2)))];
    t = t/t(end);
end

function [idx,ulocal] = mapParam(t, u)
    nt = length(t);
    nu = length(u);
    idx = zeros(1,nu);
    ulocal = zeros(1,nu);
    idx(1) = 1;
    ulocal(1) = 0;
    i = 2;
    for k=2:nu
        while u(k)>=t(i) && i<nt
            i = i+1;
        end
        i = i-1;
        idx(k) = i;
        % compute local line parameter
        ulocal(k) = (u(k)-t(i))/(t(i+1)-t(i));
    end    
end

function [Gx,Gy,Gz] = generateSurface(S, u, v)
    nu = length(u);
    nv = length(v);    
    % init outputs
    Gx = zeros(nv,nu);
    Gy = zeros(size(Gx));
    Gz = zeros(size(Gx));    
    % get polyline structs
    Gu0 = S(1);
    Gu1 = S(2);
    G0v = S(3);
    G1v = S(4);    
    % get corner points
    P00 = Gu0.p(:,1);
    P10 = Gu0.p(:,end);
    P01 = Gu1.p(:,1);
    P11 = Gu1.p(:,end);    
    % generate bilinear surface (based on Coon's bilinear interpolation)
    for j=1:nv    
        a = G0v.idx(j);
        b = G1v.idx(j);                
        P0v = G0v.p(:,a)+(G0v.p(:,a+1)-G0v.p(:,a))*G0v.h(j);
        P1v = G1v.p(:,b)+(G1v.p(:,b+1)-G1v.p(:,b))*G1v.h(j);        
        for i=1:nu            
            a = Gu0.idx(i);
            b = Gu1.idx(i);
            Pu0 = Gu0.p(:,a)+(Gu0.p(:,a+1)-Gu0.p(:,a))*Gu0.h(i);
            Pu1 = Gu1.p(:,b)+(Gu1.p(:,b+1)-Gu1.p(:,b))*Gu1.h(i);  
            % bilinear interpolation between points
            Puv = (1-v(j))*Pu0+v(j)*Pu1+(1-u(i))*P0v+u(i)*P1v - ...
                (1-u(i))*(1-v(j))*P00-v(j)*(1-u(i))*P01 - ...
                u(i)*(1-v(j))*P10-u(i)*v(j)*P11;
            Gx(j,i) = Puv(1);
            Gy(j,i) = Puv(2);
            Gz(j,i) = Puv(3);
        end
    end
end
