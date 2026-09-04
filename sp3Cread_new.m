%% ========================================================================
%  GNSS-VTEC : Ionospheric TEC Estimation from Dual-Frequency GPS Code Observations
%  ------------------------------------------------------------------------
%  Function: sp3Cread_new.m
%  Purpose : Reads a precise ephemeris SP3 (version 'a') file. Used by
%            VTEC.m for the three-day satellite orbit arcs.
%  ------------------------------------------------------------------------
%  Used by                 : Motahareh Esfandyari-Kaloukan
%  Original implementation : LaQ, MGP - TU Delft, 2003 (see credit below)
%  ========================================================================

function [sph, sp] = sp3Cread_new (file)
% SP3READ.M: read ephemerides from an SP3 file
%
% SYNTAX:
%       [sph, sp] = sp3Cread (file);
%
% DESCRIPTION:
%       file    :   input       Name of SP3 file
%       sph     :   output      Structure to hold the information from the SP3 file header
%       sp      :   output      Structure to hold the ephemeris parameters
%
% NOTE:
%       This routine can read the SP3 file with version "a" only
%
% EXAMPLE:
%       [sph, sp] = sp3Cread ('igs12003.sp3c');
%
% CREATED:
%                      -------------------
%                       28-02-2003 by LaQ
%                        MGP - TU Delft
%                      -------------------

% START:
% -------------------------------------------------------------------------

% Initialize
% fid = fopen (file);
fid = fopen (file);
if fid == -1
    error ('Cannot open SP3 file ');
end;

line = fgetl(fid);
if (~ischar (line))
    error ('Empty line in SP3 file');
end;

% Check version
if line(1:2) ~= '#a'
    error ('Not a SP3 file, version "a"');
end;

% Read header
sph.PorV = line(3);                     % Type of file (position only or with velocity)
sph.YearStart = str2num(line(4:7));     % Year of Start
sph.MonthStart = str2num(line(9:10));   % Month of Start
sph.DayStart = str2num(line(12:13));    % Day of Start
sph.HourStart = str2num(line(15:16));   % Hour of Start
sph.MinStart = str2num(line(18:19));    % Minute of Start
sph.SecStart = str2num(line(21:31));    % Second of Start
sph.NoOfEpochs = str2num(line(33:40));  % Number of Epochs
sph.DataUsed = line(41:46);             % Data used
sph.CrdSys = line(47:52);               % Coordinate system
sph.OrbitType = line(53:56);            % Orbit type
sph.Agency = line(57:60);               % Agency

line = fgetl(fid);
sph.GPSweek = str2num(line(4:7));       % GPS week
sph.SecsOfWeek = str2num(line(9:23));   % Seconds in week
sph.Interval = str2num(line(25:38));    % Interval [secs]

line = fgetl(fid);
hhh=find(line=='G');
line(hhh)=(' ');
sph.NoOfSats = str2num(line(5:6));      % Number of Satellites
%sph.SatID = zeros(1,sph.NoOfSats);
tmp = line(10:60);
for i=1:1:4
    line = fgetl(fid);
    kkk=find(line=='G');
    line(kkk)=(' ');
    tmp = [tmp line(10:60)];
end;
sph.SatID = strread(tmp,'%u',sph.NoOfSats); % List of PRNs
tmp = [];
IDMAX=max(sph.SatID);
for i=1:1:5
    line = fgetl(fid);
    tmp = [tmp line(10:60)];
end;
sph.SatAccuracy = strread(tmp,'%u',sph.NoOfSats); % Accuracy of each PRN
for i=1:1:10
    line = fgetl(fid);
end;

% Read ephemeris data
for i=1:1:sph.NoOfEpochs
    % Read time of epoch
    line = fgetl(fid);
    sp(i).Time(1) = str2num(line(4:7));
    sp(i).Time(2) = str2num(line(9:10));
    sp(i).Time(3) = str2num(line(12:13));
    sp(i).Time(4) = str2num(line(15:16));
    sp(i).Time(5) = str2num(line(18:19));
    sp(i).Time(6) = str2num(line(21:end));
    
    % Convert current epoch into GPS week and seconds in the week
    date = datenum(sp(i).Time(1), sp(i).Time(2), sp(i).Time(3)); % The serial date numbers from 1-Jan-0000 to the day of observation
    sp(i).GPSweek = fix((date - datenum(1980,1,6))/7); 
    sp(i).SecsOfWeek = (date - datenum(1980,1,6) - sp(i).GPSweek*7)*24*3600 + sp(i).Time(4)*3600 + sp(i).Time(5)*60 + sp(i).Time(6);
    
    % Read satellite postions and clock errors
    for j=1:1:sph.NoOfSats
        k=sph.SatID(j);
        line = fgetl(fid);
        ggg=find(line=='G');
        line(ggg)=(' ');
%         [ sp(i).SatID(j) sp(i).CRD(1,j) sp(i).CRD(2,j) sp(i).CRD(3,j) sp(i).Clk(j) sp(i).SD(1,j) sp(i).SD(2,j) sp(i).SD(3,j) sp(i).SClk(j)] = strread (line(2:73),' %u %f %f %f %f %f %f %f %f ');
        [ sp(i).SatID(k) sp(i).CRD(1,k) sp(i).CRD(2,k) sp(i).CRD(3,k) sp(i).Clk(k)] = strread (line(2:60),' %u %f %f %f %f ');
        % Read satellite velocities and clock rates if available
        if sph.PorV == 'V'
            line = fgetl(fid);
            [sp(i).Velo(1,k) sp(i).Velo(2,k) sp(i).Velo(3,k) sp(i).ClkRate(k)] = strread (line(5:60),'%f %f %f %f ');
        end;
    end;
end;
for i=1:1:sph.NoOfEpochs
    for j=1:IDMAX
        if sp(i).SatID(j)==0
            sp(i).CRD(:,j)=NaN;
            sp(i).Clk(j)=NaN;
            sp(i).SD(:,j)=NaN;
%             elseif sp(i).Clk(j)>=999999;
%             sp(i).CRD(:,j)=NaN;
%             sp(i).Clk(j)=NaN;
%             sp(i).SD(:,j)=NaN;
        end
    end
end
fclose(fid);

% -------------------------------------------------------------------------
% END