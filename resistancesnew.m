function [resist_out] = resistancesnew(resist_in)
    %
    %   function resistances calculates aerodynamic and boundary resistances
    %   for soil and vegetation
    %
    %   Date:       01 Feb 2008
    %   Authors:    Anne Verhoef            (a.verhoef@reading.ac.uk)
    %               Christiaan van der Tol  (tol@itc.nl)
    %               Joris Timmermans        (j_timmermans@itc.nl)
    %   Source:     Wallace and Verhoef (2000) 'Modelling interactions in
    %               mixed-plant communities: light, water and carbon dioxide', in: Bruce
    %               Marshall, Jeremy A. Roberts (ed), 'Leaf Development and Canopy Growth',
    %               Sheffield Academic Press, UK. ISBN 0849397693
    %
    %               ustar:  Tennekes, H. (1973) 'The logaritmic wind profile', J.
    %               Atmospheric Science, 30, 234-238
    %               Psih:   Paulson, C.A. (1970), The mathematical
    %               representation of wind speed and temperature in the
    %               unstable atmospheric surface layer. J. Applied Meteorol. 9,
    %               857-861
    %
    % Note: Equation numbers refer to equation numbers in Wallace and Verhoef (2000)
    %
    % Usage:
    %   [resist_out] = resistances(resist_in)
    %
    % The input and output are structures. These structures are further
    % specified in a readme file.
    %
    % Input:
    %   resist_in   aerodynamic resistance parameters and wind speed
    %
    % The strucutre resist_in contains the following elements:
    % u         =   windspeed
    % L         =   stability
    % LAI       =   Leaf Area Index

    % rbs       =   Boundary Resistance of soil                         [s m-1]
    % rss       =   Surface resistance of soil for vapour transport     [s m-1]
    % rwc       =   Within canopy Aerodynamic Resistance canopy         [s m-1]

    % z0m       =   Roughness lenght for momentum for the vegetation    [m]
    % d         =   Displacement height (Zero plane)                    [m]
    % z         =   Measurement height                                  [m]
    % h         =   Vegetation height                                   [m]

    %
    % Output:
    %   resist_out  aeorodynamic resistances
    %
    % The strucutre resist_out contains the following elements:
    % ustar     =   Friction velocity                                   [m s-1]
    % raa       =   Aerodynamic resistance above the canopy             [s m-1]
    % rawc      =   Total resistance within the canopy (canopy)         [s m-1]
    % raws      =   Total resistance within the canopy (soil)           [s m-1]

    % rai       =   Aerodynamic resistance in inertial sublayer         [s m-1]
    % rar       =   Aerodynamic resistance in roughness sublayer        [s m-1]
    % rac       =   Aerodynamic resistance in canopy layer (above z0+d) [s m-1]

    % rbc       =   Boundary layer resistance (canopy)                  [s m-1]
    % rwc       =   Aerodynamic Resistance within canopy(canopy)(Update)[s m-1]

    % rbs       =   Boundary layer resistance (soil) (Update)           [s m-1]
    % rws       =   Aerodynamic resistance within canopy(soil)          [s m-1]

    % rss       =   Surface resistance vapour transport(soil)(Update)   [s m-1]

    % uz0       =   windspeed at z0                                     [m s-1]
    % Kh        =   Diffusivity for heat                                [m2s-1]

    % load Constants
    Constants = io.define_constants();
    kappa = Constants.kappa;

    %% parameters

    Cd        =  resist_in.Cd;

    u         =  resist_in.u;
    L         =  resist_in.L;
    LAI       =  resist_in.LAI;

    % rbs       =  resist_in.rbs;
    % rss       =  resist_in.rss;
    % rwc       =  resist_in.rwc;

    z0m       =  resist_in.zo;
    d         =  resist_in.d;
    z         =  resist_in.z;
    h         =  resist_in.hc;
    w         =  resist_in.w;
    nl        =  resist_in.nl;

    p         = resist_in.p*0.1;
    % Ta        = resist_in.Ta;
    Ts        = resist_in.Ts ;
    Tsh       = Ts(1)+273.15;
    Tsu       = Ts(2)+273.15;
    po        = 101.3; % Ambient pressure (kPa)
    To        = 273.15; % Ambient Temperature (K)
    h_soil = 0.01; % roughness height of soil
    
    k_mu_soil_sh  = 1.327 * 10^(-5) * (po/p) * (Tsh/To); % KInematic viscosity 
    k_mu_soil_su  = 1.327 * 10^(-5) * (po/p) * (Tsu/To); % KInematic viscosity

    % nd = nl+1; % no of nodes for height of canopy
    dh = (h / nl) * ones(nl, 1); % incremental change of height in each layer 
    h_cum = flip([0; cumsum(dh)]); % Decreasing from top, zero at ground, hc at height of canopy 
    h_mid = (h_cum(1:end-1) + h_cum(2:end)) / 2; % At the mid of layer
    dLAI = (LAI / nl) * ones(nl, 1); % LAI in each layer.
        
    dLAI_cum = cumsum(dLAI);
    % derived parameters
    % zr: top of roughness sublayer, bottom of intertial sublayer
    zr          = 2.5 * h;                   %                            [m]
    % n: dimensionless wind extinction coefficient                       W&V Eq 33
    n           = Cd * LAI / (2 * kappa^2);      %                            []
    
    n_layer = Cd * dLAI_cum./(2*kappa^2);
    
    %% stability correction for non-neutral conditions
    % neu        = find(L >= -.001 & L <= .001);
    unst        = find(L <  -4);
    st          = find(L >  4E3);

    % stability correction functions, friction velocity and Kh=Km=Kv
    pm_z        = psim(z - d, L, unst, st);
    ph_z        = psih(z - d, L, unst, st);
    pm_h        = psim(h - d, L, unst, st);
    % ph_h       = psih(h -d,L,unst,st);
    ph_zr       = psih(zr - d, L, unst, st) .* (z >= zr) + ph_z .* (z < zr);
    phs_zr      = phstar(zr, zr, d, L, st, unst);
    phs_h       = phstar(h, zr, d, L, st, unst);
    % phih_z      = phih(z,L,unst,st);

    ustar       = max(.001, kappa * u ./ (log((z - d) / z0m) - pm_z)); %          W&V Eq 30
    

    Kh                  = kappa * ustar * (zr - d);                  %                W&V Eq 35
    % Kh                  = kappa * ustar * (zr - d) * (1/phih_z);    %W&V Eq 38
    resist_out.Kh(unst) = kappa * ustar(unst) * (zr - d) .* (1 - 16 * (h - d) ./ L(unst)).^.5; % W&V Eq 35
    resist_out.Kh(st)   = kappa * ustar(st)  * (zr - d) .* (1 + 5 * (h - d) ./ L(st)).^-1; % W&V Eq 35

    %% wind speed at height h and z0m
    uh          = max(ustar / kappa .* (log((h - d) / z0m) - pm_h), .01);
    uz          = uh * exp(n * ((h_mid / h) - 1)); 
    u_soil      = uh * exp(n * ((0.01 / h) - 1));

    Re_soil_sh     = h_soil*u_soil/k_mu_soil_sh;
    Re_soil_su     = h_soil*u_soil/k_mu_soil_su;

    St_nu_sh       = 2.46 * ((Re_soil_sh) ^ (1/4)) - log(7.4); % Stanton Number (KB-1)
    St_nu_su       = 2.46 * ((Re_soil_su) ^ (1/4)) - log(7.4); % Stanton Number (KB-1)
   
    uz0         = uh * exp(n * ((z0m + d) / h - 1));                      %       W&V Eq 32

    %% resistances

    resist_out.uz0 = uz0;
    resist_out.ustar = ustar;
    rai = (z > zr) .* (1 ./ (kappa * ustar) .* (log((z - d) / (zr - d))  - ph_z   + ph_zr)); % W&V Eq 41
    rar = 1 ./ (kappa * ustar) .*  ((zr - h) / (zr - d))      - phs_zr + phs_h; % W&V Eq 39
    % rac = calculate_rac(h,n_layer,Kh,h_mid); % W&V Eq 42
    rawc = calculate_rac_nlayer(h,n_layer,Kh,h_mid);
    raws = h * sinh(n) ./ (n * Kh) * (log((exp(n * (h) / h) - 1) / (exp(n * (h) / h) + 1)) - log((exp(n * (.01) / h) - 1) / (exp(n * (.01) / h) + 1))); % W&V Eq 43 from d of soil to height of canopy
    rbl = 70*sqrt(w./uz);  % WV Eq 31
    rbc = rbl./dLAI; % From leaf to canopy
    
    % rbc = 70 / LAI * sqrt(w ./ uz0);                        %       W&V Eq 31, but slightly different
    rbsh = St_nu_sh / (kappa * ustar); % Boundary layer resistance of soil for bare soil.
    rbsu = St_nu_su / (kappa * ustar); % Boundary layer resistance of soil for bare soil.

    resist_out.rai = rai;  
    resist_out.rar = rar;
    resist_out.rawc = rawc;  %WIthin canopy 
    resist_out.raws = raws;  % Soil (within canopy)
    resist_out.rbl = rbl;   % Boundary layer of leaf
    resist_out.rbc = rbc;   % Boundary layer of canopy integrated from leaf
    resist_out.rbsh = rbsh;   % Boundary layer of bare soil.
    resist_out.rbsu = rbsu;   % Boundary layer of bare soil.

    rac  = rai + rar + rawc;  % YT(aerodynamic , within canopy)
    rash = rai + rar + raws ; % Aerodynamic (soil, within canopy (shaded))
    rasu = rai + rar + raws;   % Boundary soil (sunlight)

    resist_out.ustar = ustar;
    resist_out.u_h = uh;
    resist_out.u_s = u_soil;
    resist_out.rac  = rac;         
    resist_out.rash = rash;        
    resist_out.rasu = rasu;   
    resist_out.L = L;

    resist_out.rac  = min(4E2, rac);         % to prevent unrealistically high resistances
    resist_out.rash = min(4E2, rash);        % to prevent unrealistically high resistances
    resist_out.rasu = min(4E2, rasu);        % to prevent unrealistically high resistances

    return

    %% subfunction pm for stability correction (eg. Paulson, 1970)
function pm = psim(z, L, unst, st)
    x           = (1 - 16 * z ./ L(unst)).^(1 / 4);
    pm          = zeros(size(L));
    pm(unst)    = 2 * log((1 + x) / 2) + log((1 + x.^2) / 2) - 2 * atan(x) + pi / 2;   %   unstable
    pm(st)      = -5 * z ./ L(st);                                      %   stable
    return

    %% subfunction ph for stability correction (eg. Paulson, 1970)
function ph = psih(z, L, unst, st)
    x           = (1 - 16 * z ./ L(unst)).^(1 / 4);
    ph          = zeros(size(L));
    ph(unst)    = 2 * log((1 + x.^2) / 2);                                %   unstable
    ph(st)      = -5 * z ./ L(st);                                      %   stable
    return

    %% subfunction ph for stability correction (eg. Paulson, 1970)
function phs = phstar(z, zR, d, L, st, unst)
    x           = (1 - 16 * z ./ L(unst)).^0.25;
    phs         = zeros(size(L));
    phs(unst)   = (z - d) / (zR - d) * (x.^2 - 1) ./ (x.^2 + 1);
    phs(st)     = -5 * z ./ L(st);
    return

function phh = phih(z,L,unst,st)
    phh         = zeros(size(L));
    phh(unst)    = (1 - 16 * z ./ L(unst)).^ (-0.5);                       % unstable 
    phh(st)      = 1 + 5 * z ./ L(st);                                      %   stable
    return

% function rac = calculate_rac(h, n, Kh, h_mid)
%     % Calculate aerodynamic resistance rac using W&V Eq. 42
%     % Inputs:
%     %   h      - canopy height (scalar)
%     %   n      - shape parameter (scalar)
%     %   Kh     - eddy diffusivity (scalar)
%     %   h_cum  - cumulative height (vector or matrix)
%     % Output:
%     %   rac    - aerodynamic resistance (same size as h_cum)
% 
%     % Compute the common terms
%     numerator = sinh(n);
%     denominator = n .* Kh;
% 
%     term1 = log((exp(n) - 1) ./ (exp(n) + 1));
% 
%     zeta = n. * (h_mid ./ h);
%     term2 = log((exp(zeta) - 1) ./ (exp(zeta) + 1));
% 
%     % Final rac expression
%     rac = h * numerator ./ denominator .* (term1 - term2);
%     return

function rac = calculate_rac_nlayer(h, n_layer, Kh, h_mid)
    % Calculate aerodynamic resistance rac using W&V Eq. 42
    % Inputs:
    %   h      - canopy height (scalar)
    %   n      - shape parameter (scalar)
    %   Kh     - eddy diffusivity (scalar)
    %   h_cum  - cumulative height (vector or matrix)
    % Output:
    %   rac    - aerodynamic resistance (same size as h_cum)

    % Compute the common terms
    numerator = sinh(n_layer);
    denominator = n_layer .* Kh;

    term1 = log((exp(n_layer) - 1) ./ (exp(n_layer) + 1));

    zeta = n_layer .* (h_mid ./ h);
    term2 = log((exp(zeta) - 1) ./ (exp(zeta) + 1));

    % Final rac expression
    rac = h * numerator ./ denominator .* (term1 - term2);
    return
