%% ========================================================================
%  GNSS-VTEC : Ionospheric TEC Estimation from Dual-Frequency GPS Code Observations
%  ------------------------------------------------------------------------
%  Script  : VTEC.m
%  Purpose : Estimates the Vertical Total Electron Content (VTEC) of the
%            ionosphere above a GPS station from dual-frequency code
%            observations. Processing chain:
%            1) Read L1/L2/P1/P2 observations (SMTReader.m).
%            2) Read three consecutive days of precise SP3 ephemerides
%               and fit cubic splines to the satellite coordinates, so
%               interpolation stays smooth across the day boundaries.
%            3) Solve the signal travel time iteratively and interpolate
%               each satellite position at the emission epoch
%               (TravelTime.m / InterpolatorFunction.m).
%            4) Compute the ionospheric pierce point of every
%               receiver-satellite ray with the single-layer model at
%               hi = 350 km.
%            5) Screen the data: elevation cut-off 3.5 deg, zero code
%               observations, low SNR.
%            6) Form the geometry-free code combination P1 - P2, correct
%               for the satellite and receiver Differential Code Biases
%               (DCB, in ns), convert the slant delay to slant TEC and
%               map it to VTEC with the cos(z') mapping function.
%            7) Plot the VTEC at the pierce points on a regional map.
%  Inputs  : ADE2003.SMT, igs15646.sp3, igs15650.sp3, igs15651.sp3
%  Requires: Mapping Toolbox
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g

spl = 299792458 ;

f1 = 1575.42e6 ; w1 = spl/f1 ;
f2 = 1227.60e6 ; w2 = spl/f2 ;

gama = f1 ^ 2 / f2 ^ 2 ;

R = 6.3781e6 ; hi = 350e3 ;

fid = fopen ('ADE2003.SMT') ;

s = SMTReader (fid) ;

[~, sp0] = sp3Cread_new ('igs15646.sp3') ;
[~, sp1] = sp3Cread_new ('igs15650.sp3') ;
[~, sp2] = sp3Cread_new ('igs15651.sp3') ;

for i = 1:size (sp0,2)
    
    xSat0 (i, :) = sp0(i).CRD(1, :) .* 1000 ;
    ySat0 (i, :) = sp0(i).CRD(2, :) .* 1000 ;
    zSat0 (i, :) = sp0(i).CRD(3, :) .* 1000 ;
    Time0 (i, :) = sp0(i).Time(:, 4:6) * [3600 60 1]' ;
    
    xSat1 (i, :) = sp1(i).CRD(1, :) .* 1000 ;
    ySat1 (i, :) = sp1(i).CRD(2, :) .* 1000 ;
    zSat1 (i, :) = sp1(i).CRD(3, :) .* 1000 ;
    Time1 (i, :) = sp1(i).Time(:, 4:6) * [3600 60 1]' ;
    
    xSat2 (i, :) = sp2(i).CRD(1, :) .* 1000 ;
    ySat2 (i, :) = sp2(i).CRD(2, :) .* 1000 ;
    zSat2 (i, :) = sp2(i).CRD(3, :) .* 1000 ;
    Time2 (i, :) = sp2(i).Time(:, 4:6) * [3600 60 1]' ;
    
end

xSat = [xSat0; xSat1; xSat2] ;
ySat = [ySat0; ySat1; ySat2] ;
zSat = [zSat0; zSat1; zSat2] ;

Time0 = Time0 - 24 * 3600 ;
Time2 = Time2 + 24 * 3600 ;

Time = [Time0; Time1; Time2] ;

for i = 1:size(xSat, 2)
    zx (:, i) = CubicSpline (Time, xSat (:, i)) ;
    zy (:, i) = CubicSpline (Time, ySat (:, i)) ;
    zz (:, i) = CubicSpline (Time, zSat (:, i)) ;
end

  xr = [-3939182.1310  3467075.3760 -3613220.8240]' ;
  
  DCBg = [01     1.998     0.259 
          02     5.090     0.015   
          03    -3.334     0.054                                          
          04    -1.977     0.068                                          
          05     0.333     0.237                                          
          06    -2.960     0.025                                          
          07     0.872     0.061                                          
          08    -3.182     0.067                                          
          09    -2.471     0.006                                          
          10    -4.346     0.109                                          
          11     1.380     0.048                                          
          12     1.742     0.042                                          
          13     1.238     0.056                                          
          14    -0.139     0.051                                          
          15     0.423     0.058                                          
          16     0.405     0.054                                          
          17     0.810     0.093                                          
          18     0.970     0.043                                          
          19     3.260     0.122                                          
          20    -1.169     0.034                                          
          21     1.658     0.042                                          
          22     5.505     0.036                                          
          23     7.151     0.115                                          
          24    -5.279     0.025                                          
          25    -1.295     0.027                                          
          26    -1.813     0.020                                          
          27    -3.469     0.043                                          
          28     0.737     0.130                                          
          29    -0.222     0.091                                          
          30    -0.739     0.095                                          
          31     2.668     0.043                                          
          32    -3.846     0.070] ; % nanosecond
  DCBr = -5.770 ; % nanosecond
  
wgs84 = referenceEllipsoid ('wgs84') ; [latr, lonr, hr] =...
    ecef2geodetic (wgs84, xr(1), xr(2), xr(3)) ; 

h = waitbar (0, 'computation of pierce point & cut off angle ...') ;

for i = 1:size(s, 2)
    t = s(i).Time(:, 4:end) * [3600 60 1]' ;
    for j = 1:s(i).NumberOfprn
        sat = s(i).PRN(j) ;
        
        [xs, tau] = TravelTime (t,xr,Time,...
            xSat(:, sat),ySat(:, sat),zSat(:, sat),...
            zx(:, sat),zy(:, sat),zz(:, sat)) ;
        
        dx = xs - xr ; u = dx ./ norm(dx) ;
        
        [E, N, U] = ecef2enu (xs(1), xs(2), xs(3), latr, lonr, hr, wgs84) ;
                
%         Vertical = atan2d  (U , ((E^2 + N^2) ^ (1/2))) ;
        Vertical = atand  (U / ((E^2 + N^2) ^ (1/2))) ;
        Zenith = 90 - Vertical ;  
        ZI = asind ((R/(R + hi)) .* sind (Zenith)) ; psy = Zenith - ZI ;
        
        S = R .* sind(psy) ./ sind(ZI) ; rPR = S .* u ;
        
        rP = xr + rPR ;
        [latP, lonP, hP] = ecef2geodetic (wgs84, rP(1), rP(2), rP(3)) ;
        
        [eP, nP, uP] = ecef2enu (rP(1), rP(2), rP(3), latr, lonr, hr, wgs84) ;
        
        IP(i).xyz(j, :) = rP' ;
        IP(i).plh(j, :) = [latP, lonP, hP] ;
        IP(i).ZI (j, :) = ZI ;
        IP(i).enu(j, :) = [eP, nP, uP] ;
        cutOffAngle(i).Angle(j, :) = Vertical ;
        
    end
    waitbar (i / size(s, 2)) 
end

close (h)

h = waitbar (0, 'deleting wrong observations') ;

for i = 1:size(s, 2)
    %%%%% Cut Off Angle < 5 degrees
    [r1, ~] = find (cutOffAngle(i).Angle < 3.5) ;
    %%%%% P1 = 0
    [r2, ~] = find (s(i).P1 <= 0) ;
    %%%%% P2 = 0
    [r3, ~] = find (s(i).P2 <= 0) ;
    %%%%% SNR
    [r4, ~] = find (s(i).SNR(:, 1) < 5) ;
    [r5, ~] = find (s(i).SNR(:, 2) < 5) ;
    [r6, ~] = find (s(i).SNR(:, 3) < 5) ;
    [r7, ~] = find (s(i).SNR(:, 4) < 5) ;
    
    s(i).P1 ([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    s(i).P2 ([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    s(i).L1 ([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    s(i).L2 ([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    
    IP(i).xyz([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    IP(i).plh([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    IP(i).ZI ([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    IP(i).enu([r1; r2; r3; r4; r5; r6; r7], :) = [] ;
    
    s(i).PRN (:, [r1; r2; r3; r4; r5; r6; r7]) = [] ;
    
    s(i).NumberOfprn = size (s(i).PRN, 2) ;
        
    waitbar (i / size(s, 2)) 
end

close (h)

h = waitbar (0, 'computation of Total Electron Content (TEC)') ;

Y = [] ;
DCBr = DCBr .* 1e-9 .* spl ;
for i = 1:size(s, 2)
    GF = s(i).P1 - s(i).P2 ;
    DCBs = DCBg (s(i).PRN, 2) .* 1e-9 .* spl ; 
    y = GF - (DCBr + DCBs) ;
    Y = [Y; y] ;
    
    waitbar (i / (2 * size (s, 2)))
    
end

N = size (Y, 1) ;

I1 = (1/(1 - gama)) .* Y ;

TEC = (I1 .* f1^2 ./ 40.3) ./ 1e16 ; % Total Electron Content

k = 1 ;
for i = 1:size(s, 2)
    for j = 1:s(i).NumberOfprn
        IP(i).VTEC(j ,:) = cosd(IP(i).ZI(j)) .* TEC(k) ;
        k = k + 1 ;
    end
    
    waitbar ((i + size (s, 2)) / (2 * size (s, 2)))
    
end

close (h)

axesm ('MapProjection', 'Robinson',...
    'Geoid', wgs84,...
    'MapLatLimit', [-50 -20],...
    'MapLonLimit', [120  150],...
    'MeridianLabel', 'on',...
    'ParallelLabel', 'on')

geoshow ('landareas.shp')

plotm (latr, lonr, 'o',...
    'MarkerFaceColor', [1 0 1],...
    'MarkerSize', 15)

h = waitbar (0, 'ploting ...') ;

for i = 1:10:size (s, 2)
    
    Time = [num2str(s(i).Time(4)), ' [h]  ',...
        num2str(s(i).Time(5)), ' [m]  ',...
        num2str(s(i).Time(6)), ' [s]  '] ;
    
    title (Time)
    
    for j = 1:size (IP(i).VTEC, 1)
   
        plotm (IP(i).plh(j, 1), IP(i).plh(j, 2), '.',...
            'MarkerFaceColor', [IP(i).VTEC(j, :)/22 0 0],...
            'MarkerEdgeColor', [IP(i).VTEC(j, :)/22 0 0])
    
    end
    
    waitbar (i / size (s, 2))
    
end

close (h)
% 
% worldmap world
% 
% geoshow ('landareas.shp')