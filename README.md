# Ionospheric VTEC from Dual-Frequency GPS Observations

MATLAB code that estimates the vertical total electron content (VTEC) of the
ionosphere above a GPS station, using dual frequency code observations and
precise satellite orbits. Written for the ADE2 station in Adelaide but works
for any station once you change the input files and the receiver coordinates.

## How it works

The geometry-free code combination P1 - P2 removes the geometry, the clocks
and the troposphere. What remains is the ionospheric delay plus the
differential code biases (DCB) of the satellite and the receiver. After
correcting for the DCBs (values in nanoseconds, taken from an IGS product),
the slant delay is converted to slant TEC and mapped to the vertical with the
single layer model, assuming a thin ionospheric shell at 350 km.

Getting the geometry right takes most of the code. Three consecutive days of
SP3 ephemerides are read and cubic splines are fitted to the satellite
coordinates, so the interpolation stays smooth at the day boundaries. For
every observation the signal travel time is solved iteratively: start with
tau = 75 ms, interpolate the satellite position at t - tau, update tau from
the geometric range, repeat until it converges. The pierce point of each ray
with the ionospheric shell is computed, observations below 3.5 degrees
elevation or with bad SNR are thrown out, and at the end the VTEC values are
drawn at their pierce points on a regional map.

The spline code is split in two on purpose. CubicSpline.m solves the
tridiagonal system once per satellite and returns the second derivatives,
and InterpolatorFunction.m evaluates the spline cheaply at any epoch. The
travel time iteration calls the interpolator thousands of times, so
refactoring the system every call would be far too slow.

## Files

- VTEC.m: main script, the whole processing chain
- SMTReader.m: reader for the SMT observation format (L1, L2, P1, P2 with LLI and SNR)
- sp3Cread_new.m: SP3-a ephemeris reader, written by LaQ at TU Delft (2003)
- CubicSpline.m: spline setup, returns the second derivatives
- InterpolatorFunction.m: spline evaluation at a query point
- TravelTime.m: iterative light-time solution
- Az_el.m: azimuth and elevation utility (the same computation is inline in VTEC.m)

## Input data

The script expects ADE2003.SMT and the ephemeris files igs15646.sp3,
igs15650.sp3 and igs15651.sp3 in the working directory. These are not in the
repository because of their size; the SP3 files can be downloaded from the
IGS/CDDIS archives. The receiver ECEF coordinates and the DCB table are set
at the top of VTEC.m and have to be adapted for other stations.

## Requirements

Mapping Toolbox (referenceEllipsoid, ecef2geodetic, ecef2enu, axesm,
geoshow, plotm).

## Notes

The elevation cut-off is set to 3.5 degrees, which is quite low; for noisy
receivers 10 degrees is safer. The color scale in the pierce point plot is
normalized by 22 TECU, so if your VTEC goes above that you should raise the
constant, otherwise MATLAB complains about color values above 1. The single
layer mapping function cos(z') gets inaccurate at very low elevations, one
more reason to be careful with the cut-off.

## Credits

sp3Cread_new.m was originally written by LaQ, MGP, TU Delft (2003), the
documentation is kept in the file header. Az_el.m follows the classical
topocent algorithm used in standard GPS toolboxes. The rest of the code is
my own work.

Motahareh Esfandyari-Kaloukan
