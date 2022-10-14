%% calculate philip's flow;

clear; close;

%% ------- figure property
set(0,'defaultAxesFontSize',15);
set(0,'defaultAxesFontName','Arial')
set(0,'defaultlegendFontName','Arial')
set(0,'defaulttextinterpreter','latex');

L = 1.0; % normalize parameters, L can be fixed 

% if you plot the coutour of the velocity field, plotf == 1
plotf = 1;
% if you want to evaluate the value on the meniscus, calculate 
% ZpU = exp(1i*(0:0.01:1)*pi); % x_men = zzeta(ZpU); w_men = wzeta(ZpU); 
% You can plot like plot(x_men,w_men);




%% ----- Modify here. You can use appropriate value delta, and q
delta = 0.3354070926037155953025603594142;
q = 0.1102644608366191436576997375596;
mesh_size = 0.0125;



%% define prime function
dv = [delta*1i -delta*1i]; qv = [q q];
thinf = dv - qv.^2./(conj(dv));

w = skprod_vec(dv,qv,5);
chizeta = @(zeta) w(zeta,thinf(1))./w(zeta,thinf(2));
zzeta = @(zeta) -1i*L/pi*log(chizeta(zeta));
B = real(zzeta(-1)); H = imag(zzeta(qv(1)*1i+dv(1)));

%% solve the boundary value problem
bv0 = @(zeta) -H*real(zzeta(zeta));
bv1 = @(zeta) 0; bv2 = @(zeta) 0;
bvfun = {bv0,bv1,bv2};

dth = pi/400;
Hzeta = @(zeta) Schwarz_integ_multiply_onlyC0(zeta,bvfun,dv,qv,dth);
wzeta = @(zeta) imag(zzeta(zeta)).*(H-imag(zzeta(zeta))/2) + imag(Hzeta(zeta));


if plotf == 1
[X,Y] = meshgrid([-1:mesh_size:1]);
mask = ((X.^2+Y.^2)<1+1e-6)+((X-real(dv(1))).^2+(Y-imag(dv(1))).^2>(qv(1)-1e-6).^2) + ...
        ((X-real(dv(2))).^2+(Y-imag(dv(2))).^2>(qv(2)-1e-6).^2) == 3;
dthp = pi/30; Zpp = exp(1i*(0+dthp/2:dthp:2*pi-dthp/2)');    
Z = X(mask) + 1i*Y(mask);    
[X,Y] = meshgrid([-1:0.0025:1]);
mask = ((X.^2+Y.^2)<1+1e-6) + ((X.^2+Y.^2)>0.99) + (Y>-0.01) == 3;
Z = [Z;X(mask) + 1i*Y(mask);];

brh = 1i*transpose(0:0.02:imag(dv(1))-qv(1));
ep = 1e-5;
Z = [Z;Zpp;qv(1)*Zpp+dv(1);transpose(-1:0.01:1);brh-ep;brh+ep];
% add
Z = [Z;qv(1)*Zpp*1.01+dv(1);transpose(-1:0.01:1)+1i*0.01];
zcd = zzeta(Z);
lenz = numel(zcd);
nn = 100;
klz = floor(lenz/nn);
wzcd = zeros(size(Z));
disp("Calculating the flow...");
for ll = 1:nn
    disp([ll+"/"+nn]);
    wzcd((ll-1)*klz+1:ll*klz) = wzeta(transpose(Z((ll-1)*klz+1:ll*klz)));
end
wzcd(nn*klz+1:end) = wzeta(transpose(Z(nn*klz+1:end)));

dm = 0.0025;
[xq,yq] = meshgrid(0:dm:L, 0:dm:H);
vq = griddata(real(zcd),imag(zcd),wzcd,xq,yq);

%% plot ---

figure()
colorb = [0 0.4470 0.7410];
contour(xq,yq,vq,50,'LineColor',colorb);
hold on;
contour(-xq,yq,vq,50,'LineColor',colorb);
contour(xq,2*H-yq,vq,50,'LineColor',colorb);
contour(-xq,2*H-yq,vq,50,'LineColor',colorb);
contour(2*L-xq,yq,vq,50,'LineColor',colorb);
contour(2*L-xq,2*H-yq,vq,50,'LineColor',colorb);
contour(-2*L+xq,yq,vq,50,'LineColor',colorb);
contour(-2*L+xq,2*H-yq,vq,50,'LineColor',colorb);

plot(-2*L+B:0.01:-B,0*(-2*L+B:0.01:-B),'k','Linewidth',2.2);
plot((L+B:-0.01:B),0*(L+B:-0.01:B),'k','Linewidth',2.2);
plot(-2*L+B:0.01:-B,0*(-2*L+B:0.01:-B)+2*H,'k','Linewidth',2.2);
plot((L+B:-0.01:B),0*(L+B:-0.01:B)+2*H,'k','Linewidth',2.2);
plot(-B:0.01:B,0*(-B:0.01:B),'b--','Linewidth',2.2);
plot(-B:0.01:B,0*(-B:0.01:B)+2*H,'b--','Linewidth',2.2);
plot(L+(B:0.01:L),0*(B:0.01:L),'b--','Linewidth',2.2);
plot(L+(B:0.01:L),0*(B:0.01:L)+2*H,'b--','Linewidth',2.2);
plot(-2*L+(0:0.01:B),0*(0:0.01:B),'b--','Linewidth',2.2);
plot(-2*L+(0:0.01:B),0*(0:0.01:B)+2*H,'b--','Linewidth',2.2);

title("Contour plot $$w(x,y)$$",'Interpreter','latex','Fontsize',20)
axis equal;
xlim([-2*L-0.1,2*L+0.1]);
ylim([-0.1,2*H+0.1]);

end




function F = Schwarz_integ_multiply_onlyC0(z,bvfun,dv,qv,dth)
% vector representation of Villat formula

M = length(qv);
% ---------------------------------------------
skv = skprod_vec(dv, qv, 5);
angz = angle(z);
% ---- alternate trapezoidal rule: --- 
labels = zeros(M+1,numel(z));
labels(1,:) = abs(abs(z)-1) < 1e-5;
addp = bvfun{1}(z).*labels(1,:);
for j = 1:M
    labels(j+1,:) = abs(abs(z-dv(j)) - qv(j)) < 1e-5;
    addp = addp + bvfun{j+1}(z).*labels(j+1,:);
    angz = angz.*(1-labels(j+1,:)) + angle((z-dv(j))/qv(j)).*labels(j+1,:);
end

% -- find angle and add value

th = transpose(0:dth:2*pi-dth) + angz;
Zp = exp(1i*th);
th0 = transpose(0-dth/2:dth:2*pi-dth/2) + angz;
Zp0 = exp(1i*th0);


%% check
bv{1} = bvfun{1}(Zp);
for j = 1:M
    bv{j+1} = bvfun{j+1}(qv(j)*Zp+dv(j));
end

kp0 = log(skv(z,Zp0).^2./Zp0);
Kerp{1} = diff(real(kp0) + 1i*unwrap(imag(kp0)));
%for j = 1:M
%    cz = 1./conj(z); cz(cz==Inf) = 1e+8;
%    kp1 = log(skv(z,qv(j)*Zp0 + dv(j))) + log(conj(skv(cz,qv(j)*Zp0 + dv(j))));
%    Kerp{j+1} = diff(real(kp1) + 1i*unwrap(imag(kp1)));
%end

F = 1/(2*pi*1i)*sum(bv{1}.*Kerp{1});

%for j = 1:M
%    F = F - 1/(2*pi*1i)*sum(bv{j+1}.*Kerp{j+1});
%end

F = F + addp;


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



