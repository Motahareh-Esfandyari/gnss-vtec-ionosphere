%% ========================================================================
%  GNSS-VTEC : Ionospheric TEC Estimation from Dual-Frequency GPS Code Observations
%  ------------------------------------------------------------------------
%  Function: SMTReader.m
%  Purpose : Reads a GPS observation file in the SMT format (17 header
%            lines, then per epoch: time, epoch flag, number of PRNs and
%            list of PRNs, followed by one line per satellite with
%            L1, L2, P1, P2 and their LLI and SNR flags). Returns a
%            structure array with one element per epoch.
%  Inputs  : fid - file identifier from fopen()
%  Outputs : s   - structure array (Time, EpochFlag, NumberOfprn, PRN,
%            L1, L2, P1, P2, LLI, SNR)
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [s] = SMTReader (fid)


txt = textscan ( fid , '%s' , 'delimiter' , '\n' , 'HeaderLines' , 17 ) ;

h = waitbar (0, 'reading the observation file ...') ;

k = 1 ; o = 1 ; i = 1 ;

A = str2num (txt{1, 1}{1, 1}) ;
        
s(1).Time        = A (1:6) ;
s(1).EpochFlag   = A (7) ;
s(1).NumberOfprn = A (8) ;
s(1).PRN         = A (9:end) ;
while i ~= size(txt{1, 1}, 1)
    for sat = 1:s(k).NumberOfprn
        i = i + 1 ;
        %%% L1
        l1 = [] ;
        for j = 1:size(txt{1, 1}{i, 1}, 2)
            if strcmp(txt{1, 1}{i, 1}(1, j), '.') ~= 1
                l1 = [l1,txt{1, 1}{i, 1}(1, j)] ;
            else
                j1 = j ;
                break
            end
        end
        l1 = [l1,txt{1, 1}{i, 1}(1, j1:j1 + 3)] ; 
        s(k).L1(sat, :) = str2num(l1) ;
        %%% LLI
        lli = str2num(txt{1, 1}{i, 1}(1, j1 + 4)) ;
        if isempty(lli) == 1
            lli = NaN ;
        end
        s(k).LLI(sat, 1) = lli ;
        %%% SNR
        snr = str2num(txt{1, 1}{i, 1}(1, j1 + 5)) ;
        if isempty(snr) == 1
            snr = NaN ;
        end
        s(k).SNR(sat, 1) = snr ;
        %%% L2
        l2 = [] ;
        for j = j1 + 6:size(txt{1, 1}{i, 1}, 2)
            if strcmp(txt{1, 1}{i, 1}(1, j), '.') ~= 1
                l2 = [l2,txt{1, 1}{i, 1}(1, j)] ;
            else
                j2 = j ;
                break
            end
        end
        l2 = [l2,txt{1, 1}{i, 1}(1, j2:j2 + 3)] ;
        s(k).L2(sat, :) = str2num(l2) ;
        %%% LLI
        lli = str2num(txt{1, 1}{i, 1}(1, j2 + 4)) ;
        if isempty(lli) == 1
            lli = NaN ;
        end
        s(k).LLI(sat, 2) = lli ;
        %%% SNR
        snr = str2num(txt{1, 1}{i, 1}(1, j2 + 5)) ;
        if isempty(snr) == 1
            snr = NaN ;
        end
        s(k).SNR(sat, 2) = snr ;
        %%% P1
        p1 = [] ;
        for j = j2 + 6:size(txt{1, 1}{i, 1}, 2)
            if strcmp(txt{1, 1}{i, 1}(1, j), '.') ~= 1
                p1 = [p1,txt{1, 1}{i, 1}(1, j)] ;
            else
                j3 = j ;
                break
            end
        end
        p1 = [p1,txt{1, 1}{i, 1}(1, j3:j3 + 3)] ;
        s(k).P1(sat, :) = str2num(p1) ;
        %%% LLI
        lli = str2num(txt{1, 1}{i, 1}(1, j3 + 4)) ;
        if isempty(lli) == 1
            lli = NaN ;
        end
        s(k).LLI(sat, 3) = lli ;
        %%% SNR
        snr = str2num(txt{1, 1}{i, 1}(1, j3 + 5)) ;
        if isempty(snr) == 1
            snr = NaN ;
        end
        s(k).SNR(sat, 3) = snr ;
        %%% P2
        p2 = [] ;
        for j = j3 + 6:size(txt{1, 1}{i, 1}, 2)
            if strcmp(txt{1, 1}{i, 1}(1, j), '.') ~= 1
                p2 = [p2,txt{1, 1}{i, 1}(1, j)] ;
            else
                j4 = j ;
                break
            end
        end
        p2 = [p2,txt{1, 1}{i, 1}(1, j4:j4 + 3)] ;
        s(k).P2(sat, :) = str2num(p2) ;
        %%% LLI
        lli = str2num(txt{1, 1}{i, 1}(1, j4 + 4)) ;
        if isempty(lli) == 1
            lli = NaN ;
        end
        s(k).LLI(sat, 4) = lli ;
        %%% SNR
        snr = str2num(txt{1, 1}{i, 1}(1, j4 + 5)) ;
        if isempty(snr) == 1
            snr = NaN ;
        end
        s(k).SNR(sat, 4) = snr ;
        
        waitbar (i / size(txt{1, 1}, 1))
    end
    if i == size(txt{1, 1}, 1)
        break
    end
    k = k + 1 ; i = i + 1 ;
    A = str2num (txt{1, 1}{i, 1}) ;
        
    s(k).Time        = A (1:6) ;
    s(k).EpochFlag   = A (7) ;
    s(k).NumberOfprn = A (8) ;
    s(k).PRN         = A (9:end) ;
    
    waitbar (i / size(txt{1, 1}, 1))
    
end

close (h)