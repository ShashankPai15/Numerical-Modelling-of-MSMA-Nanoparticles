function Ps = clausius_clapeyron(T, P0, Hlv, Rgas, Tb)
% CLAUSIUS_CLAPEYRON  Saturation vapor pressure of the workpiece metal.
%
% Dimensionally-consistent Clausius-Clapeyron relation (corrected form of
% the paper's eq.13 -- see Modelling Note 1 at the top of
% microEDM_nanoparticle_model.m).
%
%   Ps = P0 * exp( Hlv*(T-Tb) / (Rgas*T*Tb) )

    exponent = Hlv*(T - Tb) / (Rgas*T*Tb);
    exponent = min(exponent, 700);   % overflow guard (exp(700) ~ 1e304)
    Ps = P0*exp(exponent);
end
