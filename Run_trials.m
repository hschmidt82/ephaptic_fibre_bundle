% This script generates the main results (Figure 8) shown in Schmidt, Hahn,
% Deco, Knoesche; PLOS Computational Biology (2020).

close all; 
clear all;

Np = 5; % number of trials for same parameter set
Nm = 4; % number of different stimulus durations
No = 10; % number of different stinulus intensities
Nq = 4; % number of different bundle radii

Ntot = Nm*Np*No*Nq; % total number of trials

% arrays to record trial outputs:
delay = NaN(Ntot,1);
delay2 = delay;
delay0 = delay;
delay02 = delay;
lat = delay;
lat0 = delay;

% map parfor loop variable t onto trial variables m,o,p,q:
mvec = repmat([1:Nm],[1 Np*No*Nq]);
vp = ones(Nm,1)*[1:Np];
pvec = repmat(vp(:)',[1 No*Nq]);
vo = ones(Nm*Np,1)*[1:No];
ovec = repmat(vo(:)',[1 Nq]);
vq = ones(Nm*Np*No,1)*[1:Nq];
qvec = vq(:)';

% define parameters to be passed on to trial function:
for t = 1:Ntot
    
    m = mvec(t);
    p = pvec(t);
    o = ovec(t);
    q = qvec(t);
    
    pars(t).dens = o/No;
    pars(t).Deltat = m;
    pars(t).pref = 1;
    pars(t).radbund = 4*q/Nq;
    pars(t).gr = 0.8;
    pars(t).rho = 0.8;
    pars(t).sigrat = 3/(1-pars(t).rho);
    
    pars0(t) = pars(t);
    pars0(t).pref = 0;
    
end

parfor t = 1:Ntot
    
    m = mvec(t);
    p = pvec(t);
    o = ovec(t);
    q = qvec(t);
    
    % trial with ephaptic coupling:
    temp = fSpikerun(pars(t));
    delay(t) = temp(2);
    delay2(t) = temp(3);
    lat(t) = temp(1);
    
    % trial without ephaptic coupling:
    temp = fSpikerun(pars0(t));
    delay0(t) = temp(2);
    delay02(t) = temp(3);
    lat0(t) = temp(1);
    
end

for t = 1:Ntot
    m = mvec(t);
    p = pvec(t);
    o = ovec(t);
    q = qvec(t);
    dmat(m,p,o,q) = delay(t);
    dmat2(m,p,o,q) = delay2(t);
    dmat0(m,p,o,q) = delay0(t);
    dmat02(m,p,o,q) = delay02(t);
    lmat(m,p,o,q) = lat(t);
    lmat0(m,p,o,q) = lat0(t);
    
end

save('Results.mat','dmat','dmat0','dmat2','dmat02','lmat','lmat0')

% END OF SCRIPT