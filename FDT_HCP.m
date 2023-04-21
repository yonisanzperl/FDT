clear all;
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Tenet/TENET/'];
addpath(genpath(path3));

%% Version with sliding windows...

tasks={'REST1';'EMOTION';'GAMBLING';'WM';'LANGUAGE';'MOTOR';'RELATIONAL';'SOCIAL'};

N=62;
NSUB=970;
Tau=3;
indexN=[1:31 50:80];
sigma=0.01;

Isubdiag = find(tril(ones(N),-1));

% Parameters of the data
TR=0.72;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

load empirical_dbs80;
load SC_dbs80HARDIFULL.mat;
C = SC_dbs80HARDI;
C=C(indexN,indexN);
C = C/max(max(C));

Cefftask2=zeros(NSUB,N,N);
Cefftask=zeros(size(tasks,1),N,N);
FDTdeviation=zeros(size(tasks,1),NSUB);
permap=zeros(size(tasks,1),NSUB,N);
Onsagerdiff=zeros(size(tasks,1),NSUB);

for xx=1:size(tasks,1)
    load(['hcp1003_' tasks{xx} '_LR_dbs80.mat']);
    nsub=1;
    for sub=1:1003
        if isstruct(subject{sub})
            if size(subject{sub}.dbs80ts,2)>175 % one subject has only in EMOTION
                subject{nsub}=subject{sub};
                nsub=nsub+1;
            end
        end
    end
    
    for sub=1:NSUB
        fprintf('%d \t %d \n',xx,sub)
        ts2=subject{sub}.dbs80ts;
        ts=ts2(indexN,1:176);
        clear signal_filt;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            %             signal_filt(seed,:)=(filtfilt(bfilt,afilt,ts(seed,:)));
        end
        ts2=ts(:,10:end-10);
        Tm=size(ts2,2);
        FCemp=corrcoef(ts2');
        COVemp=cov(ts2');
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
            [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff,sigma);
            COVtausim=expm((Tau*TR)*A)*COVsimtotal;
            COVtausim=COVtausim(1:N,1:N);
            %         COVtausim=COVtausim/max(max(COVtausim));
            for i=1:N
                for j=1:N
                    sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
                end
            end
            COVtausim=COVtausim.*sigratiosim;
            errorFC(iter)=mean(mean((FCemp-FCsim).^2));
            errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));
            
            if mod(iter,100)<0.1
                errornow=mean(mean((FCemp-FCsim).^2));                
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
                        Cnew(i,j)=Cnew(i,j)+0.0002*(FCemp(i,j)-FCsim(i,j)) ...
                            +0.0002*(COVtauemp(i,j)-COVtausim(i,j));
                        if Cnew(i,j)<0
                            Cnew(i,j)=0;
                        end
                    end
                end
            end
            Cnew = Cnew/max(max(Cnew));
        end
        Ceff=Cnew;
        Cefftask2(sub,:,:)=Ceff;
        
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff,sigma);
        invA=inv(A);
        for i=1:2*N
            for j=1:2*N
                hh=zeros(2*N,1);
                hh(j)=1;
                xepsilon=-invA*hh;
                chi(i,j)=abs((2*COVsimtotal(i,j)/sigma^2)-xepsilon(i));
                chi2(i,j)=abs(xepsilon(i));
            end
        end
        chij=mean(chi(1:N,1:N))./mean(chi2(1:N,1:N));
        FDTdeviation(xx,sub)=mean(chij);
        permap(xx,sub,:)=chij(1:N);
        % Onsager
        Onsagerdiff(xx,sub)=mean(mean(abs(Ceff-Ceff')));
    end
    Cefftask(xx,:,:)=squeeze(mean(Cefftask2));
end

figure(1);
boxplot([FDTdeviation(1,:)' FDTdeviation(2,:)' FDTdeviation(3,:)' FDTdeviation(4,:)' FDTdeviation(5,:)' FDTdeviation(6,:)' FDTdeviation(7,:)' FDTdeviation(8,:)']);
for xx=2:8
    a=FDTdeviation(1,:);
    b=FDTdeviation(xx,:);
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],50000,0.01,'ttest');
    pFDT(xx)=stats.pvals(2);
end

figure(2);
boxplot([Onsagerdiff(1,:)' Onsagerdiff(2,:)' Onsagerdiff(3,:)' Onsagerdiff(4,:)' Onsagerdiff(5,:)' Onsagerdiff(6,:)' Onsagerdiff(7,:)' Onsagerdiff(8,:)']);
for xx=2:8
    a=Onsagerdiff(1,:);
    b=Onsagerdiff(xx,:);
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],50000,0.01,'ttest');
    pAsym(xx)=stats.pvals(2);
end

save results_HCP.mat Cefftask Onsagerdiff permap FDTdeviation;