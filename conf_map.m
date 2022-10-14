
%% plot conformal mapping
% After "determine_params.m", we can use this code to check the map.
% 


clear; close;
%% figure property
set(0,'defaultAxesFontSize',15);
set(0,'defaultAxesFontName','Arial')
set(0,'defaultlegendFontName','Arial')
set(0,'defaulttextinterpreter','latex');

L = 1.0; % normalize parameters, L can be fixed 

%% -----
delta = 0.3354070926037155953025603594142;
q = 0.1102644608366191436576997375596;

dv = [delta*1i -delta*1i]; qv = [q q];
thinf = dv - qv.^2./(conj(dv));

%% define prime function
w = skprod_vec(dv,qv,5);
chizeta = @(zeta) w(zeta,thinf(1))./w(zeta,thinf(2));
zzeta = @(zeta) -1i*L./pi * log(chizeta(zeta));

B = zzeta(-1); H = imag(zzeta(qv(1)*1i+dv(1)));
[X,Y] = meshgrid([-1:0.1:1]);
mask = ((X.^2+Y.^2)<1)+((X-real(dv(1))).^2+(Y-imag(dv(1))).^2>qv(1).^2) + ...
        ((X-real(dv(2))).^2+(Y-imag(dv(2))).^2>qv(2).^2) == 3;
Zp = exp(1i*pi*(0:0.01:2)); 
ZpU = exp(1i*pi*(0:0.01:1)); 


ep = 1.0e-3;
db = -imag(dv(1))+qv(1):0.01:imag(dv(1))-qv(1);

%% --- plot conformal mapping
% plot triply connected domain
figure()
subplot(1,2,1);
plot(Zp,'k-','LineWidth',1.5); 
hold on;
fill(real([ZpU,-1:0.1:1]),imag([ZpU,-1:0.1:1]),[0.9,0.9,0.9])
fill(real(qv(1)*Zp+dv(1)),imag(qv(1)*Zp+dv(1)),[1,1,1])
plot(qv(1)*Zp+dv(1),'b-','LineWidth',1.5);
plot(qv(2)*Zp+dv(2),'r-','LineWidth',1.5);
plot((-1:0.01:1),0*(-1:0.01:1),'m-','LineWidth',1.5);
plot(X(mask)+1i*Y(mask),'k.','MarkerSize',0.4);
plot(ep+db*0,db,'k--','LineWidth',1.5);
plot(-ep+db*0,db,'c--','LineWidth',1.5);
title("$$\zeta$$-plane",'Fontsize',20,'Interpreter','latex');
axis square;

% plot physical domain
subplot(1,2,2);
fill([-L,L,L,-L],[0,0,H,H],[0.9,0.9,0.9])
hold on;
plot(zzeta(Zp),'k-','LineWidth',1.5); 
plot(zzeta(qv(1)*Zp+dv(1)),'b-','LineWidth',1.5);
plot(zzeta(qv(2)*Zp+dv(2)),'r-','LineWidth',1.5);
plot(zzeta((-1:0.01:1)),'m.','LineWidth',1.5);
plot(zzeta(X(mask)+1i*Y(mask)),'k.','MarkerSize',0.04);
plot(zzeta(ep+db*0+1i*db),'k--','LineWidth',1.5);
plot(zzeta(-ep+db*0+1i*db),'c--','LineWidth',1.5);
title("$$\mathcal{Z}$$-plane",'Fontsize',20,'Interpreter','latex');
axis equal;
xlim([-L-0.1,L+0.1]);
ylim([-H-0.1,H+0.1]);




function cf = cf_delta_qv(delta1,q1,a,b,H)
dv = [delta1 -delta1];
qv = [q1 q1];
thinf = dv - qv.^2./(conj(dv));
zeta0 = delta1 + q1*1i;
%% define prime function
w = skprod_vec(dv,qv,6);
chizeta = @(zeta) w(zeta,thinf(1))./w(zeta,thinf(2));
zzeta = @(zeta) -1i*a/pi*log(chizeta(zeta));
cf(1) = abs(zzeta(-1)-b);
cf(2) = abs(zzeta(zeta0)-1i*H);
cf = sum(cf);
end





%%%% -------------------------------- Schottoky-Klein prime function
function wf = skprod_vec(dv, qv, L)
%SKPROD is the product formula for the S-K prime function.
%
%   wf = skprod(dv, qv);
%   w = wf(z, alpha);
% 
% This function computes a truncated half Schottky group necessary
% to calculate the truncated product formula for the Schottky-Klein
% prime function. It returns a function handle which evaluates the prime
% function.
%
% -----
% Input:
%   dv = a vector of circle centers.
%   qv = a vector of circle radii.
%   L = (optional) truncation level of the product formula (default L=4).
%
% Output:
%   wf = a function handle to the prime function with signature
%     w = wf(z, alpha),
%   where z is an array of points at which to evaluate the function, and
%   alpha is a scalar parameter value.
%
% The "hat" version of the product formula, where the zero and pole have
% been factored out, may be accessed via
%    w_hat = wf(z, alpha, 'hat').
%
% -------
% Example:
% 
%   >> dv = [0.5, 0.5i];
%   >> qv = [0.1, 0.1];
%   >> wf = skprod(dv, qv, 6);
%   >> w = wf(-0.5-0.5i, 1);
%   >> w^2
%
%   ans =
%   
%              2.39754812001042 +      1.76164377385124i
%
% This value may be checked against
%   D. G. Crowdy and J. S. Marshall, "Computing the Schottky-Klein prime
%   function on the Schottky double of planar domains," CMFT 7 (2007) no.
%   1, 293-308.
%
% The prime function is documented in
%   H. Baker, Abelian Functions and the Allied Theory of Theta Functions,
%   Cambridge University Press, Cambridge, 1897, 1995.

% Copyright Everett Kropf, 2015
% 
% skprod is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% skprod is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with skprod.  If not, see <http://www.gnu.org/licenses/>.

if nargin < 3
    lmax = 4;
else
    lmax = L;
end

if numel(dv) ~= numel(qv)
    error('It is expected the numbers of centers and radii be equal.')
end

%Group generators.
m = numel(dv);
ngen = 2*m;
gens = cell(ngen, 1);
for j = 1:m
    th = [qv(j)^2 - abs(dv(j))^2, dv(j); -conj(dv(j)), 1]/qv(j);
    gens{j} = th;
    gens{j+m} = [th(4), -th(3); -th(2), th(1)];
end

% Group setup.
nhg = (ngen*(ngen - 1).^(0:lmax-1))/2;
grp = cell(sum(nhg), 1);
glvl = [0, cumsum(nhg(1:end-1))];

% Permutation-y matrix
perm = zeros(ngen, ngen-1);
for n = 1:ngen
    perm(n,:) = find(1:ngen ~= mod(n + ngen/2 - 1, ngen) + 1);
end

% Prep for search.
word = ones(1, lmax);
nptr = ones(1, lmax);
node = cell(lmax, 1);
lastg = zeros(1, lmax);

recur = @(f,varargin) f(f, varargin{:});
cmp = @(x,y) recur(@(f,x,y) ...
    x(1) < y(1) || x(1) == y(1) && f(f, x(2:end), y(2:end)), x, y);

function nfill(word)
    wlvl = numel(word);
    if wlvl > 1 && isempty(node{wlvl-1})
        nfill(word(1:end-1));
    end
    node{wlvl} = node{wlvl-1}*gens{word(wlvl)};
end

% Depth first search.
lvl = 0;
while true
    % Next deeper.
    lvl = lvl + 1;
    if lvl > 1
        word(lvl) = perm(word(lvl - 1), nptr(lvl));
        % Is inverse already in group?
        if lastg(lvl) > 0
            winv = mod(word(lvl:-1:1) + ngen/2 - 1, ngen) + 1;
            wchk = cmp(word, winv);
        else
            wchk = true;
        end
        if wchk
            lastg(lvl) = lastg(lvl) + 1;
            if isempty(node{lvl})
                nfill(word(1:lvl));
            end
            grp{glvl(lvl)+lastg(lvl)} = node{lvl};
        end
    else
        word(1) = nptr(1);
        node{1} = gens{word(1)};
        if word(1) <= ngen/2
            lastg(1) = lastg(1) + 1;
            grp{lastg(1)} = node{1};
        end
    end
    if lvl < lmax
        continue
    end
    
    while lvl > 0
        % Go back!
        node{lvl} = [];
        lvl = lvl - 1;
        
        % Try to turn.
        if nptr(lvl+1) < ngen - (lvl > 0)
            % Can turn.
            nptr(lvl+1) = nptr(lvl+1) + 1;
            break
        else
            % Can't turn.
            nptr(lvl+1) = 1;
        end
    end
    
    if (lvl == 0 && nptr(1) == 1) || all(lastg == nhg)
        break
    end
end

% Evaluate.
function w = skeval(z, alpha, hat)
    %if numel(alpha) ~= 1
    %    error('The parameter alpha must be a scalar.')
    %end
    
    if nargin > 2 && strcmp(hat, 'hat')
        dohat = true;
    else
        dohat = false;
    end
    
    if isa(z, 'double')
        islarge = abs(z) > 2^1000;
    elseif isa(z, 'single')
        islarge = abs(z) > 2^120;
    else
        error('Argument z must be single or double precision.')
    end
    atinf = isinf(z);
    
    if ~isinf(alpha(1)) && ~dohat
        w = z - alpha;
    else
        w = complex(ones(size(z)));
    end
    
    function fprod(th)
        thjz = (th(1)*z + th(3))./(th(2)*z + th(4));
        if ~isinf(alpha)
            thja = (th(1)*alpha + th(3))./(th(2)*alpha + th(4));
            w = w.*(thjz - alpha).*(thja - z)...
                ./(thjz - z)./(thja - alpha);
        else
            w = w.*(th(1)/th(2) - z)./(thjz - z);
        end
    end
    cellfun(@fprod, grp)
    
    if ~isinf(alpha(1))
        w(atinf) = inf;
        % Normalize failed large number input.
        w(islarge & isnan(w)) = inf;
    else
        % No pole.
        w(atinf) = 1;
        w(islarge & isnan(w)) = 1;
    end
end

wf = @skeval;

end



