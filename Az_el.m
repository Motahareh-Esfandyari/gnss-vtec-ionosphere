%% ========================================================================
%  GNSS-VTEC : Ionospheric TEC Estimation from Dual-Frequency GPS Code Observations
%  ------------------------------------------------------------------------
%  Function: Az_el.m
%  Purpose : Computes the azimuth, elevation angle and geometric range
%            from a receiver to one or more satellites. Transforms the
%            ECEF baseline into the local topocentric ENU frame at the
%            receiver and derives the angles from the ENU components.
%            Utility function; the same computation is done inline in
%            VTEC.m via ecef2enu.
%  Inputs  : XP - receiver ECEF position [m]; XS - satellite ECEF
%            positions [m], one column per satellite
%  Outputs : Az [deg, 0-360], El [deg], D - range [m]
%  Note    : Follows the classical topocent algorithm used in standard
%            GPS toolboxes; adapted to use ecef2geodetic.
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [Az, El, D] = Az_el(XP,XS);

         
dx=[ XS(1,:)-XP(1); XS(2,:)-XP(2); XS(3,:)-XP(3) ];

D = sqrt(dx(1,:).^2+dx(2,:).^2+dx(3,:).^2);
dtr = pi/180;
% [phi,lambda,h] = XYZ2FIL(XP(1),XP(2),XP(3),'WGS84');
[phi,lambda,h] = ecef2geodetic (XP(1), XP(2), XP(3), referenceEllipsoid('wgs84')) ;

phi = phi .* (pi/180) ; lambda = lambda .* (pi/180) ;

cl = cos(lambda); sl = sin(lambda);
cp = cos(phi); sp = sin(phi);
R = [-sl -sp*cl cp*cl;
      cl -sp*sl cp*sl;
       0    cp   sp];
Topov = R'*dx;
E = Topov(1,:);
N = Topov(2,:);
U = Topov(3,:);
h_dis = sqrt(E.^2+N.^2);

for iii=1:length(E)

if h_dis(iii) < 1.e-20
   Az(iii) = 0;
   El(iii) = 90;
else
   Az(iii) = atan2(E(iii),N(iii))/dtr;
   El(iii) = atan2(U(iii),h_dis(iii))/dtr;
end
if Az(iii) < 0
   Az(iii) = Az(iii)+360;
end
end

