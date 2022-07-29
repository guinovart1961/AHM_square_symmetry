function [LB , UB]=Hashin79(z1,z2,z3,z4,lam)

Eam=z1(1);
Etm=z1(2);
Gam=z1(3);
Gtm=z1(4);
Pam=z1(5);
Ptm=z1(6);

Eaf=z2(1);
Etf=z2(2);
Gaf=z2(3);
Gtf=z2(4);
Paf=z2(5);
Ptf=z2(6);

c11m=z3(1);
c12m=z3(2);
c13m=z3(3);
c33m=z3(4);
c44m=z3(5);
c66m=(c11m-c12m)/2;
km=(c11m+c12m)/2;
                    
c11f=z4(1);
c12f=z4(2);
c13f=z4(3);
c33f=z4(4);
c44f=z4(5);
c66f=(c11f-c12f)/2;
kf=(c11f+c12f)/2;


%*****************************************************************************************
% Cotas de Hashin 1979
%******** Axiales
%Modulo Axial de Young 
LEa=Eam*(1-lam)+Eaf*lam+ 4*(Pam-Paf)^2*(1-lam)*lam/...
   (2*(1-lam)/(c11f+c12f)+(2*lam/(c11m+c12m))+(1/Gtm));

UEa=Eam*(1-lam)+Eaf*lam+ 4*(Pam-Paf)^2*(1-lam)*lam/...
   (2*(1-lam)/(c11f+c12f)+(2*lam/(c11m+c12m))+(1/Gtf));

%Poisson Axial
LPa=Pam*(1-lam)+Paf*lam+(Paf-Pam)*...
   (2/(c11m+c12m)-2/(c11f+c12f))*(1-lam)*lam/...
   (2*(1-lam)/(c11f+c12f)+(2*lam/(c11m+c12m))+...
   (1/Gtf));

UPa=Pam*(1-lam)+Paf*lam+(Paf-Pam)*...
(2/(c11m+c12m)-2/(c11f+c12f))*(1-lam)*lam/...
(2*(1-lam)/(c11f+c12f)+(2*lam/(c11m+c12m))+ ...
(1/Gtm));

%Modulo  Axial de cizalladura
LGa=Gam+lam/(1/(Gaf-Gam)+(1-lam)/(2*Gam));

UGa=Gaf+(1-lam)/(1/(Gam-Gaf)+lam/(2*Gaf));

%Modulo de Bulk
Lk=km+ lam/(1/(kf-km)+(1-lam)/(km+Gtm));
Uk=kf+ (1-lam)/(1/(km-kf)+(lam)/(kf+Gtf));

%************* Transversales
%****LGt

%beta1=km/(km+2*Gtm);
beta1=1/(3-4*Pam);
beta2=kf/(kf+2*Gtf);
gamma=Gtf/Gtm;
alfa=(beta1-gamma*beta2)/(1+gamma*beta2);
ro=(gamma+beta1)/(gamma-1);
LGt=Gtm+lam/(1/(Gtf-Gtm)+(km+2*Gtm)*(1-lam)/...
      (2*Gtm*(km+Gtm)));

UGt=Gtf+(1-lam)/(1/(Gtm-Gtf)+(kf+2*Gtf)*(lam)/...
                      (2*Gtf*(kf+Gtf)));
      
%UGt=Gtm*(1+((1+beta1)*lam)/(ro-lam*(1+(3*beta1^2*(1-...
 %   lam)^2)/(alfa*(lam)^3+1))));
Lm=1+(4*Lk*UPa^2)/LEa;
Um=1+(4*Uk*LPa^2)/UEa;
%LEt=(4*Lk*LGt)/(Lk+Lm*LGt);

LEt=4/(1/LGt+1/Lk+4*UPa^2/LEa);
%UEt=4/(1/UGt+1/Lk+4*UPa^2/LEa);

UEt=(4*Uk*UGt)/(Uk+Um*UGt);
LPt=(Lk-Lm*UGt)/(Lk+Lm*UGt);
UPt=(Uk-Um*LGt)/(Uk+Um*LGt);


%LG1=1/(Gtf-Gtm);
%LG2=(1-lam)*(km+2*Gtm)/(2*Gtm*(km+Gtm));
%LGt=Gtm+lam/(LG1+LG2);
%****UGt
%gm1=(Gtf/Gtm);
%beta1=km/(km+2*Gtm);
%beta2=kf/(kf+2*Gtf);
%alf1=(beta1-gm1*beta2)/(1+gm1*beta2);
%ro1=(gm1+beta1)/(gm1-1);
%de1=3*beta1^2*(1-lam)^2/(alf1*(lam)^3+1);

%UGt=Gtm*( 1+ (1+beta1)*lam/( ro1 - lam*(1+de1) ));
%********
%Lm=1+4*Lk*UPa/LEa;
%Um=1+4*Uk*UPa/UEa;

%******** Et
%LEt=4*Lk*LGt/(Lk+Lm*LGt);
%UEt=4*Uk*UGt/(Uk+Um*UGt);

%*************Pt
%LPt=(Lk-Lm*UGt)/(Lk+Lm*UGt);
%UPt=(Uk-Um*LGt)/(Uk+Um*LGt);

% DATOS EXPERIMENTALES DE ABOUDI PAG. 60
UB=[UEa UEt UGa UGt UPa UPt Uk];

LB=[LEa LEt LGa LGt LPa LPt  Lk];

