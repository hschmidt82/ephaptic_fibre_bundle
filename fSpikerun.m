% This function numerically integrates the spike propagation model for a
% given set of parameters.
function [output] = spikerun(pars)

% convert parameters passed on to function:
dens = pars.dens;
Deltat = pars.Deltat;
pref = pars.pref;
radbund = pars.radbund;
sigrat = pars.sigrat;
gr = pars.gr;
rho = pars.rho;

% parameters of bundle geometry:
Nf = 1e4; % number of fibres / axons
L = 1e2; % length of fibres in mm;
dx = 1e-1; % spatial discretisation in mm;
x = [0:dx:L]; % discretised space;

% sample diameter distribution:
xx = [0:1e-3:5];
a = 0.2; b = 0.25;
gf = @(x) max((1/b^2).*(x-a).*exp((a-x)./b),0);
gmax = max(gf(xx)); xmax = max(xx);
for m = 1:Nf;
    miss = 1;
    while miss
        rnd1 = xmax.*rand;
        rnd2 = gmax.*rand;
        if rnd2<gf(rnd1)
            miss = 0;
            df(m) = rnd1;
        end
    end
end

% convert fibre diameters to velocity distribution:
alpha = 5;
vf = alpha.*df;

% parameter to scale the effect of ephaptic field on velocity:
Vfit = 180;

% ephaptic field:

Vmax = 100; % maximum amplitude of action potential (from resting potential)

x1 = 0.3; % time to peak of action potential
x3 = 2; % duration of action potential
x2 = x3-x1;

% coupling function (for piecewise linear approximation of action potential):
eph = @(t,ddist,vel) pref.*(ddist./sum(df.^2)).*(sigrat*gr^2*rho/2).*...
    ((Vmax./(x1.*vel)).*(sqrt(radbund.^2 + (t-x).^2) - sqrt(radbund.^2 + (t-x-x1.*vel).^2) - abs(t-x) + abs(t-x-x1.*vel))-...
    (Vmax./(x2.*vel)).*(sqrt(radbund.^2 + (t-x-x1.*vel).^2) - sqrt(radbund.^2 + (t-x-x3.*vel).^2) - abs(t-x-x1.*vel) + abs(t-x-x3.*vel)));

% simulation parameters:
T = 1e3; % simulation time in ms
dt = 1e-2; % time step in ms
T = round(T/dt);
T2 = T/10;

% set up spike initiation times for each fibre:
ST = zeros(Nf,1); % spike initiation time
[~,ind] = sort(rand(Nf,1)); % shuffled fibre indices
Nfrac = round(dens*Nf); % fraction of active fibres
for m = 1:Nfrac
    ST(ind(m),1) = ceil((Deltat/dt)*rand)+1; % assign random time to each active fibre
end
STA = ST;

% spike position:
pos = zeros(Nf,1);
vmat = vf';
vmat2 = vmat;
mask = pos; mask2 = pos;

% Jansen-Rit variables:
y0 = zeros(1,T+1);
y1 = y0; y2 = y0; y3 = y0; y4 = y0; y5 = y0;
vq = -3.65.*ones(T+1,1); rq = 0.09.*ones(T+1,1); syn = zeros(T+1,1);

% Jansen-Rit parameters:
A = 3.25; B = 22;
a = 100; b = 50;
e0 = 2.5; r = 0.56;
v0 = 6; C1 = 135;
C2 = 0.8*C1; C3 = 0.25*C1; C4 = C3;

% numerical integration routine:
t = 1; % start at first time step
while (t<T)&&(min(STA(ind)-ST(ind))==0)
    
    % identify currently active fibres:
    mask = 1.*((pos>0)&(pos<=L));
    
    % identify position, diameter, and velocity of currently active fibres:
    spikepos = pos(mask==1);
    spikedf = df(mask==1);
    spikevmat = vmat2(mask==1);
    
    % compute EP along fibre bundle:
    feph2 = 0;
    for m = 1:sum(mask)
        feph2 = feph2 + eph(spikepos(m),spikedf(m)^2,spikevmat(m))';
    end
    
    % compute EP at position of each action potential:
    indi = floor(spikepos./dx);
    indl = spikepos./dx - indi;
    feph = (1-indl).*feph2(indi+1) + indl.*feph2(indi+2);
    
    % compute velocity of spikes:
    spikev = vmat(mask==1)./(1 +1.*(1/Vfit).*feph);
    
    % compute new position of spikes:
    pos(mask==1) = pos(mask==1) + dt.*spikev;
    
    % initiate spikes whose spike initiation times have been reached:
    pos(ST==t) = dx;
    
    % record arrival times of spikes that have reached the distal end (pos>L)
    STA(mask2(:)-mask(:)>0)=t;
    
    % compute number of spikes that have reached the distal end (for Jansen-Rit model)
    finst = sum(1.*(mask2(:)-mask(:)>0));
    
    % reset mask array for next time step:
    mask2 = mask;
    
    % compute velocity of centre of mass of spikes (required to compute EP)
    vmat2(mask==0) = vmat(mask==0);
    vmat2(mask==1) = (1-dt).*vmat2(mask==1) + dt.*spikev;

    % response of Jansen-Rit model to spike volleys:
    y0(t+1) = y0(t) + 1e-3*dt.*y3(t);
    y3(t+1) = y3(t) + 1e-3*dt.*(A*a*2*e0/(1+exp(r*(v0-(y1(t)-y2(t))))) - 2*a*y3(t) - a^2*y0(t));
    y1(t+1) = y1(t) + 1e-3*dt.*y4(t);
    y4(t+1) = y4(t) + 1e-3*dt.*(A*a*((1e6/(Nf*dt))*finst+C2*2*e0/(1+exp(r*(v0-C1*y0(t))))) - 2*a*y4(t) - a^2*y1(t));
    y2(t+1) = y2(t) + 1e-3*dt.*y5(t);
    y5(t+1) = y5(t) + 1e-3*dt.*(B*b*C4*2*e0/(1+exp(r*(v0-C3*y0(t)))) - 2*b*y5(t) - b^2*y2(t));
    
    t = t+1; % advance to next time step
end

% record maximum response time of Jansen-Rit model as proxy for latency:
[~,tmax] = max(y1-y2);

% compute mean and STD of delay distribution:
STA = STA(find(STA>0));
ST = ST(find(ST>0));
delay = mean(STA(:)-ST(:)); % mean
delay2 = std(STA(:)-ST(:)); % STD

output = [tmax,delay,delay2];

% END OF FUNCTION