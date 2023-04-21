clear all
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Tenet/TenetFCtau/'];
addpath(genpath(path3));

N=62; %% numero de regiones 
NSUB=18; %% numero de sujetos

indexN=[1:31 50:80];  %% Cortical areas
Tau=1;
sigma=0.01;

epsFC=0.0002;
epsFCtau=0.00004;
maxC=0.1;

load empirical_sleep_W_DK.mat;
load laufs_sleep.mat;

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

C = sc(indexN,indexN);  %% anatomie..DTI tractography
C = C/max(max(C))*maxC;

%%
% Wake
for nsub=1:NSUB
    nsub
    ts=TS_W{nsub};  % fMRI wakefulness
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCW(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=C;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,A]=Lineal_int(Cnew,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsim;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.01
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N  %% learning
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceff=Cnew;
    CeffW(nsub,:,:)=Ceff;
    [FCsim,COVsim,A]=Lineal_int(Cnew,sigma);
    fittFC_W(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsim;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_W(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end

% N3
load empirical_sleep_N3_DK.mat;

for nsub=1:NSUB
    nsub
    ts=TS_N3{nsub};
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCN3(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=C;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,A]=Lineal_int(Cnew,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsim;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.01
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceff=Cnew;
    CeffN3(nsub,:,:)=Ceff;
    [FCsim,COVsim,A]=Lineal_int(Cnew,sigma);
    fittFC_N3(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsim;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_N3(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end

%% FDT violation
for nsub=1:NSUB
    Ceff=squeeze(CeffW(nsub,:,:));
    [FCsim,COVsim,A]=Lineal_int(Ceff,sigma);
    invA=inv(A);
    for i=1:N
        for j=1:N
            hh=zeros(N,1);
            hh(j)=1;
            xepsilon=-invA*hh;
            chi(i,j)=abs((2*COVsim(i,j)/sigma^2)-xepsilon(i));
            chi2(i,j)=abs(xepsilon(i));
        end
    end
    chij=(chi./chi2);
    FDTW(nsub,:,:)=chij;
    FDTWm(nsub)=mean(chij(:));
    perFDTW(nsub,:)=mean(chij);
    perCeffW(nsub,:)=mean(Ceff);
    perFCW(nsub,:)=mean(squeeze(FCW(nsub,:,:)));
end

for nsub=1:NSUB
    Ceff=squeeze(CeffN3(nsub,:,:));
    [FCsim,COVsim,A]=Lineal_int(Ceff,sigma);
    invA=inv(A);
    for i=1:N
        for j=1:N
            hh=zeros(N,1);
            hh(j)=1;
            xepsilon=-invA*hh;
            chi(i,j)=abs((2*COVsim(i,j)/sigma^2)-xepsilon(i));
            chi2(i,j)=abs(xepsilon(i));
        end
    end
    chij=(chi./chi2);
    FDTN3(nsub,:,:)=chij;
    FDTN3m(nsub)=mean(chij(:));       
    perFDTN3(nsub,:)=mean(chij);
    perCeffN3(nsub,:)=mean(Ceff);
    perFCN3(nsub,:)=mean(squeeze(FCN3(nsub,:,:)));
end

figure(1)
boxplot([FDTWm' FDTN3m']);
a=FDTWm;
b=FDTN3m;
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],50000,0.01,'ttest');
stats.pvals

figure(2);
boxplot([fittFC_W' fittCVtau_W' fittFC_N3' fittCVtau_N3']);

save results_Lineal_FDT_sleepDK.mat perFCW perFCN3 perFDTW perFDTN3 perCeffW perCeffN3 FDTWm FDTN3m fittFC_N3 fittFC_W fittCVtau_N3 fittCVtau_W FDTW FDTN3 FCW FCN3;
