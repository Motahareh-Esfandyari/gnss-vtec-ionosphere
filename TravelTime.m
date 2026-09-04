%% ========================================================================
%  GNSS-VTEC : Ionospheric TEC Estimation from Dual-Frequency GPS Code Observations
%  ------------------------------------------------------------------------
%  Function: TravelTime.m
%  Purpose : Iterative light-time solution. Starting from tau = 75 ms,
%            interpolates the satellite position at the emission epoch
%            t - tau with the precomputed cubic splines and updates
%            tau = |xs - xr| / c until convergence (1e-10 s).
%  Inputs  : t - reception time [s]; xr - receiver ECEF position [m];
%            Time, xSat, ySat, zSat - satellite ephemeris samples;
%            zx, zy, zz - spline second derivatives (CubicSpline.m)
%  Outputs : xs - satellite ECEF position at emission [m]; tau - travel
%            time [s]
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [xs,tau] = TravelTime (t,xr,Time,xSat,ySat,zSat,zx,zy,zz)



c = 299792458 ;

tau0 = 0.075 ;
error = 1 ;
while error > 1e-10
    
    [xInt] = InterpolatorFunction (Time, xSat, t - tau0, zx) ;
    [yInt] = InterpolatorFunction (Time, ySat, t - tau0, zy) ;
    [zInt] = InterpolatorFunction (Time, zSat, t - tau0, zz) ;
    
    xs = [xInt yInt zInt]' ;
    
    tau = (1/c) * norm (xs - xr) ;
    
    error = abs (tau - tau0) ; tau0 = tau ;
    
end

end