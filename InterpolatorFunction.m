%% ========================================================================
%  GNSS-VTEC : Ionospheric TEC Estimation from Dual-Frequency GPS Code Observations
%  ------------------------------------------------------------------------
%  Function: InterpolatorFunction.m
%  Purpose : Evaluates the natural cubic spline at the query point x0,
%            given the sample points (x, y) and the second derivatives z
%            precomputed by CubicSpline.m.
%  Inputs  : x, y - sample points ; x0 - query point ; z - from CubicSpline
%  Outputs : yInt - interpolated value at x0
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [yInt] = InterpolatorFunction (x, y, x0, z)



[r1,~] = find (x == x0) ; h  = diff (x) ;

if isempty(r1) ~= 1
    
    yInt = y (r1,:) ;
    
else
    
    [r2,~] = find (x < x0)  ; i = r2 (end) ;

    yInt = z(i + 1)/(6 * h(i)) * (x0 - x(i))^3 + z(i)/(6 * h(i)) * ...
        (x(i + 1) - x0)^3 + (y(i + 1)/h(i) - (h(i) * z(i + 1))/6) * ...
        (x0 - x(i)) + (y(i)/h(i) - (h(i) * z(i))/6) * (x(i + 1) - x0) ;
    
end