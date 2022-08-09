function EE=ingenieriles3_Hashin79(z)

c11=z(1);
c12=z(2);
c13=z(3);
c33=z(4);
c44=z(5);
c66=z(6);
k=(c11+c12)/2;
mp=z(7);

Pa=c13/(2*k); % Poisson axial efectivo
Ea=c33-4*k*(Pa)^2; % Ea efectivo
Ga=c44;
Gt=mp;
m=1+4*k*(Pa)^2/Ea;
Et=4*k*Gt/(k+m*Gt); % Et efectivo
Pt=(k-m*Gt)/(k+m*Gt); % Poisson t efectivo

%******************************************************
EE=[Ea Et Ga Gt Pa Pt];
