function [resist_out] = resistancesHF(resist_in)
%==========================================================================
% resistancesHF - Compute aerodynamic and boundary layer resistances
%                 for high-frequency land-atmosphere modeling using Harman
%                 and Finigan , 2007 and 2008
%
% INPUT:
%   resist_in - struct containing input parameters:
%       Cd   - drag coefficient
%       u    - wind speed (m/s)
%       L    - Obukhov length (m)
%       LAI  - leaf area index
%       z    - measurement height (m)
%       h    - canopy height (m)
%       w    - leaf width (m)
%       nl   - number of canopy layers
%       p    - pressure (kPa)
%       Ta   - air temperature (K)
%
% OUTPUT:
%   resist_out - struct containing resistance components and ustar
%==========================================================================

% Load constants
Constants = io.define_constants();
kappa = Constants.kappa;

%--------------------------------------------------------------------------
% Unpack input
Cd = resist_in.Cd;
u  = resist_in.u;
L  = resist_in.L;
LAI = resist_in.LAI;
z  = resist_in.z;
hc = resist_in.hc;
w  = resist_in.w;
nl = resist_in.nl;
p  = resist_in.p;
Ta = resist_in.Ta;
Ts = resist_in.Ts;
Tsh = 273.15 + Ts(1); %shaded
Tsu = 273.15 + Ts(2); %sunlit

% Constants and empirical coefficients
c2 = 0.5;
beta_neutral_max = 0.35;
po = 101.3;       % Ambient pressure (kPa)
To = 273.15;      % Ambient temperature (K)
z0ms = 0.01;      % Roughness length of soil (m)
D0 = 2.5e-5;  % molecular diffusivity of vapor in air (m2/s) at 298 K (Standard)

% Kinematic viscosity of air (should be at soil using soil temperature)
kmu_sh = 1.327e-5 * (po / p) * (Tsh / To) ^ (1.81);  %shaded
kmu_su = 1.327e-5 * (po / p) * (Tsu / To) ^ (1.81);  %sunlit
%--------------------------------------------------------------------------
% Load lookup tables for ψ_hat (roughness sublayer corrections)
lookup_psihat = "psihat.nc";
dtLgridM = ncread(lookup_psihat, 'dtLgridM');
zdtgridM = ncread(lookup_psihat, 'zdtgridM');
psigridM = ncread(lookup_psihat, 'psigridM');
dtLgridH = ncread(lookup_psihat, 'dtLgridH');
zdtgridH = ncread(lookup_psihat, 'zdtgridH');
psigridH = ncread(lookup_psihat, 'psigridH');

%--------------------------------------------------------------------------
% Vertical canopy structure
nd = nl + 1;
dh = (hc / nl) * ones(nl, 1);
h_cum = flip([0; cumsum(dh)]);             % Height from top to ground
hmid = (h_cum(1:end-1) + h_cum(2:end)) / 2;
dLAI = (LAI / nl) * ones(nl, 1);
dLAI_cum = cumsum(dLAI);

%--------------------------------------------------------------------------
% Calculate canopy length scale
Lc = 4 * hc / LAI;

% Compute neutral beta and limit it
c_beta = kappa^2 / log((hc + z0ms) / z0ms)^(-2);
beta_neutral = min(sqrt(c_beta + 0.3 * LAI), beta_neutral_max);

% Avoid near-zero L
if abs(L) <= 0.1
    L = 0.1;
end

LcL = Lc / L;

% Compute beta = u*/u(h)
beta = GetBeta(beta_neutral, LcL);

% Compute Prandtl and Schmidt numbers
[Pr, Sc] = GetPrSc(beta_neutral, beta_neutral_max, LcL);

%--------------------------------------------------------------------------
% Displacement height (d) from top of canopy
dt = beta^2 * Lc * (1 - exp(-0.25 * LAI / beta^2));
dt = min(hc, dt);
d = hc - dt;

% Limit Obukhov length to avoid extreme cases, the phim, phic, psim, psic
% are only valid from -2 to 1, Foken, 2006. 
zeta = (z - d) / L;
if zeta >= 0
    zeta = min(1.0, max(zeta, 0.01));
else
    zeta = max(-2.0, min(zeta, -0.01));
end
L = (z - d) / zeta;

%--------------------------------------------------------------------------
% Stability functions
phim = phi_m_monin_obukhov((hc - d) / L);
phic = phi_c_monin_obukhov((hc - d) / L);

% Roughness sublayer correction coefficients
c1_m = (1 - kappa / (2 * beta * phim)) * exp(0.5 * c2); %For momentum
c1_c = (1 - (Sc * kappa) / (2 * beta * phic)) * exp(0.5 * c2); %For heat,vapor

% Deviations from log-law due to roughness sublayer
dv_m = - psi_m_monin_obukhov((z - d)/(L)) + ...
         psi_m_monin_obukhov((hc - d)/L) + ...
         c1_m * (LookupPsihat((z - hc)/(hc - d), (hc - d)/L, zdtgridM, dtLgridM, psigridM) - ...
         LookupPsihat((hc - hc) / (hc - d), (hc - d) / L, zdtgridM, dtLgridM, psigridM)) + kappa/beta;

% dv_c = - psi_c_monin_obukhov((z - d)/(hc - d)) + ...
%          psi_c_monin_obukhov((hc - d)/L) + ...
%          c1_c * (LookupPsihat((z - hc)/(hc - d), (hc - d)/L, zdtgridH, dtLgridH, psigridH) - ...
%          LookupPsihat((hc - hc) / (hc- d), (hc- d) / L, zdtgridH, dtLgridH, psigridH));

%--------------------------------------------------------------------------
% Wind profile and turbulence characteristics
zlog = log((z - d) / (hc - d));
ustar = u * kappa / (zlog + dv_m);
% ustar_su = u * kappa / (log((z) / (z0ms)) - psi_m_monin_obukhov((z)/(L)) + ...
%          psi_m_monin_obukhov((z0ms)/L)); %No canopy effect 
uh = ustar / beta;
lm = 2 * beta^3 * Lc; %Mixing length, formula (6, H&F 2007)
uz = uh * exp((hmid - hc) / (lm / beta)); % Wind speed at each layer
Kch = lm * ustar / Sc; % Diffusivity at canopy top. 
% usu = ustar_su/kappa * (log((0.02) / (z0ms)) - psi_m_monin_obukhov((z)/(L)) + ...
%          psi_m_monin_obukhov((z0ms)/L)); % very close to surface
us = uh * exp((0.02 - hc) / (lm / beta)); % %Wind speed at soil height. 

%--------------------------------------------------------------------------
% Aerodynamic and boundary layer resistances
rar = (1 / (kappa * ustar)) * (zlog + dv_m); % Resitance from the top of the canopy to the flux tower height.
rawc = (1 / (beta * ustar)) * ...
       (exp(-(hmid - hc) / (lm / beta)) - exp(-(hc- hc) / (lm / beta))); % Within canopy resistance, second term is zero, because to canopy height.
raws = (1 / (beta * ustar)) * ...
       (exp(-(z0ms - hc) / (lm / beta)) - exp(-(hc- hc) / (lm / beta)));% Within canopy resistance for the soil part.
% rasu = ((1/kappa*ustar)) * (log((z) / (z0ms)) - psi_m_monin_obukhov((z)/(L)) + ...
%          psi_m_monin_obukhov((z0ms)/L)); % Aerodynamic resistance from bare soil to reference height

rar  = min(400, rar);
rawc = min(400, rawc);
raws = min(400, raws);

uz   = max(uz,0.01);
us  = max(us,0.01);


% Boundary layer part
rbl = 70 * sqrt(w ./ uz);   % Leaf boundary layer resistance
% rbc = rbl ./ dLAI;          % Leaf to canopy
rbc = rbl;

% Soil boundary layer
Re_su = z0ms * us / kmu_su; %Reynold's number
Re_sh = z0ms * us / kmu_sh;

% For now I cancel it , because if windspeed too low. 
St_nu_su     = 2.46 * ((Re_su) ^ (1/4)) - log(7.4);  %Stanton Number, sunlit
St_nu_sh     = 2.46 * ((Re_sh) ^ (1/4)) - log(7.4);  %Stanton Number, shaded
rbsu         = St_nu_su / (kappa * ustar);          
rbsh         = St_nu_sh / (kappa * ustar);

deltasu = 5* z0ms /(sqrt(Re_su)); %Just to check
deltash = 5* z0ms /(sqrt(Re_sh));
% 
% rbsu = deltasu/D0; %Flat plate style
% rbsh = deltash/D0; %Flat plate style

%--------------------------------------------------------------------------
% Final resistances (I think this should be including bl)
rac = rar+rawc;
rash = rar+raws; 
rasu = rar+raws; % Same for sunlit and sun shaded.

 
rbsu = max(4, min(400, rbsu));
rbsh = max(4, min(400, rbsh));

%--------------------------------------------------------------------------
% Pack output
resist_out.ustar = ustar;
% resist_out.ustarsu = ustar_su;
resist_out.u_s = us;
% resist_out.ush = ush;
resist_out.u_h = uh;
resist_out.d = d;
resist_out.rar = rar;
resist_out.rawc = rawc;
resist_out.raws = raws;
resist_out.rasu = rasu;
resist_out.rbl = rbl;
resist_out.rbc = rbc;
resist_out.rbsu = rbsu;
resist_out.rbsh = rbsh;
resist_out.rac = rac;
resist_out.rash = rash;
resist_out.rasu = rasu;
resist_out.deltasu = deltasu;
resist_out.deltash = deltash;
resist_out.lm = lm;
resist_out.beta = beta;
resist_out.kmu_sh = kmu_sh;
resist_out.kmu_su = kmu_su;
resist_out.L = L;
end

