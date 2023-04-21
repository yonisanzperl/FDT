clear all;
load results_Lineal_FDT_sleepDK.mat;

%% group change of Hierarchy
N=62;

pwm=mean(perFDTW);
psm=mean(perFDTN3);
[sopwm indpwm]=sort(pwm-psm,'descend');

figure(1);
sopwm=mean(perFDTW(:,indpwm));
sopsm=mean(perFDTN3(:,indpwm));
sopws=std(perFDTW(:,indpwm));
sopss=std(perFDTN3(:,indpwm));
shadedErrorBar(1:N,sopwm,sopws/2/sqrt(18),'b',0.7)
hold on;
shadedErrorBar(1:N,sopsm,sopss/2/sqrt(18),'r',0.7)
plot(sopwm-psm(indpwm),'k','LineWidth',2);


pwm=mean(perCeffW);
psm=mean(perCeffN3);
[sopwm indpwm]=sort(pwm-psm,'descend');
figure(2)
soWceffm=mean(perCeffW(:,indpwm));
soN3ceffm=mean(perCeffN3(:,indpwm));
soWceffs=std(perCeffW(:,indpwm));
soN3ceffs=std(perCeffN3(:,indpwm));
shadedErrorBar(1:N,sopwm,sopws/2/sqrt(18),'b',0.8)
hold on;
shadedErrorBar(1:N,soWceffm,soWceffs/2/sqrt(18),'c',0.4)
shadedErrorBar(1:N,soN3ceffm,soN3ceffs/2/sqrt(18),'r',0.4)

pwm=mean(perFCW);
psm=mean(perFCN3);
[sopwm indpwm]=sort(pwm-psm,'descend');
figure(3)
soWfcm=mean(perFCW(:,indpwm));
soN3fcm=mean(perFCN3(:,indpwm));
soWfcs=std(perFCW(:,indpwm));
soN3fcs=std(perFCN3(:,indpwm));
shadedErrorBar(1:N,sopwm,sopws/2/sqrt(18),'b',0.8)
hold on;
shadedErrorBar(1:N,soWfcm,soWfcs/2/sqrt(18),'c',0.4)
shadedErrorBar(1:N,soN3fcm,soN3fcs/2/sqrt(18),'r',0.4)

%% Top in hierarchy change

%% Top Analysis based on group...
pwm=mean(perFDTW);
psm=mean(perFDTN3);
[sopwm indpwm]=sort(pwm-psm,'descend');
pww=perFDTW(:,indpwm);
pss=perFDTN3(:,indpwm);
ppw=[];
pps=[];
for i=1:62
    pw=mean(pww(:,1:i),2);
    ps=mean(pss(:,1:i),2);
    ppw=[ppw pw];
    pps=[pps ps];
    p(i)=ranksum(pw,ps);
end
figure(4)
boxplot([ppw(:,1) pps(:,1) ppw(:,5) pps(:,5) ppw(:,10) pps(:,10) ]);

%% Bottom analysis
ppw=[];
pps=[];
for i=1:61
    pw=mean(pww(:,end-i:end),2);
    ps=mean(pss(:,end-i:end),2);
    ppw=[ppw pw];
    pps=[pps ps];
    p(i)=ranksum(pw,ps);
end
figure(5)
boxplot([ppw(:,1) pps(:,1) ppw(:,20) pps(:,20) ]);

%%

