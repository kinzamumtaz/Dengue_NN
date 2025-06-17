function dy = Dengue3( t,y);

dy=zeros(17,1);


piH=26; 
piV=2600;
eta1=1.4;
eta2=2.8; 
muH=(1/65)/365; 
muV=26/365;
phi1=73/365;
phi2=73/365;
gamma1=52/365;
gamma2=52/365;
deltaH1=0.18/365;
deltaH2=0.18/365;
deltaV=0/365;
omega1=12/365;
omega2=12/365;
u1=0.6;
u2=0.3;
alphaH1=20/365;
alphaH2=20/365;
alphaV1=1000/365;   
alphaV2=800/365;

NH=y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14);
NV=y(15)+y(16)+y(17);

lambdaH1 = alphaH1*(eta1*y(13)+y(3))/NH;
lambdaV1 = alphaV1*(y(16))/NH;

lambdaH2 = alphaH2*(eta2*y(7)+y(9))/NH;
lambdaV2 = alphaV2*(y(17))/NH;

dy(1)=piH-lambdaV1*y(1)-lambdaV2*y(1)-muH*y(1);

dy(2)=lambdaV1*y(1)-(phi1 + muH)*y(2);

dy(3)=phi1*y(2)-(gamma1 + muH+ deltaH1)*y(3);

dy(4)=gamma1*y(3)-(muH+omega1)*y(4);

dy(5)=omega1*y(4)-(u1*lambdaV2+(1-u1)+muH)*y(5);

dy(6)=u1*lambdaV2*y(5)-(phi2+muH)*y(6);

dy(7)=phi2*y(6)-(gamma2 + muH+ deltaH2)*y(7);

dy(8)=lambdaV2*y(1)-(phi2 + muH)*y(8);

dy(9)=phi2*y(8)-(gamma2 + muH+ deltaH2)*y(9);

dy(10)=gamma2*y(9)-(muH+omega2)*y(10);

dy(11)=omega2*y(10)-(u2*lambdaV1+(1-u2)+muH)*y(11);

dy(12)=u2*lambdaV1*y(11)-(phi1+muH)*y(12);

dy(13)=phi1*y(12)-(gamma1 + muH+ deltaH1)*y(13);

dy(14)=gamma1*y(13) +gamma2*y(7) +(1-u1)*y(5) +(1-u2)*y(11) -muH*y(14);

dy(15)=piV-lambdaH1*y(15)-lambdaH2*y(15)-muH*y(15);

dy(16)=lambdaH1*y(15)-(deltaV+muV)*y(16);

dy(17)=lambdaH2*y(15)-(deltaV+muV)*y(17);
end

