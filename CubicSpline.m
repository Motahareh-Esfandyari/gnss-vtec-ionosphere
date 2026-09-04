%% ========================================================================
%  GNSS-VTEC : Ionospheric TEC Estimation from Dual-Frequency GPS Code Observations
%  ------------------------------------------------------------------------
%  Function: CubicSpline.m
%  Purpose : Natural cubic spline setup. Solves the tridiagonal system
%            once and returns the vector of second derivatives z, which
%            InterpolatorFunction.m then uses to evaluate the spline at
%            any query point. Splitting setup and evaluation avoids
%            refactoring the system for every epoch, which matters here
%            because the travel-time iteration calls the interpolator
%            thousands of times.
%  Inputs  : x, y - sample points (sorted by x)
%  Outputs : z    - second derivatives at the nodes
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [z] = CubicSpline (x,y)



h  = diff (x) ;

dy = diff (y) ;

b  = dy./h ;

for i = 2:size (h,1)
    
    U (i - 1 , :) = 2 .* (h(i - 1 , :) + h(i , :)) ;
    V (i - 1 , :) = 6 .* (b(i , :) - b(i - 1 , :)) ;
    
end

n = size (U,1) ;

B = zeros (n,n) ; B (1,1:2) = [U(1) h(2)] ; B (n,n - 1:n) = [h(n) U(n)] ;

for i = 2:n - 1
    
    B (i , i - 1:i + 1) = [h(i) U(i) h(i + 1)] ;
    
end

z = inv(B) * V ;

z = [0 ; z ; 0] ;

end