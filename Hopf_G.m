clear all
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Tenet/TenetFCtau/'];
addpath(genpath(path3));

N=90;
Tau=1;
sigma=0.01;

load DataSleepW_N3.mat;
load empirical_sleep_N3.mat;

Isubdiag = find(tril(ones(N),-1));

% Parameters of the data
TR=2.08;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

C = SC;
C = C/max(max(C))*0.2;


%% Wake
NSUB=15;
for nsub=1:NSUB
    nsub
    ts=TS_N3{nsub};
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        %         signal_filt(seed,:)=(filtfilt(bfilt,afilt,ts(seed,:)));
    end
    ts2=ts(:,10:end-10);
    Tm=size(ts2,2);
    FCemp2(nsub,:,:)=corrcoef(ts2');
end
FCemp=squeeze(mean(FCemp2));
Cnew=C;
iter2=1;
GG=0.:0.01:5;
for g=GG
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(g*Cnew,f_diff,sigma);
    error(iter2)=mean(mean((FCemp-FCsim).^2));
    iter2=iter2+1;
end
[aux ind]=min(error);
g=GG(ind)