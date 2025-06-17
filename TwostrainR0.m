PiH=26; 
PiV=2600;
eta1=1.4;
eta2=28; 
muH=(1/65)/365; 
muV=26/365;
phi1=73/365;
phi2=73/365;
gamma1=52/365;
gamma2=52/365;
delta1=0.18/365;
delta2=0.18/365;
deltaV=0/365;
omega1=12/365;
omega2=12/365;
u1=0.6;
u2=0.3;
alphaH1=10/365;
alphaH2=10/365;
alphaV1=1000/365;   
alphaV2=800/365;



R_0(1)=sqrt((phi1*alphaH1*muH*PiV*alphaV1)/((phi1+muH)*(gamma1+delta1+muH)*(deltaV+muV)*PiH*muV))
R_0(2)=sqrt((phi2*alphaH2*muH*PiV*alphaV2)/((phi2+muH)*(gamma2+delta2+muH)*(deltaV+muV)*PiH*muV))
