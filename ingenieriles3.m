function EE=ingenieriles3(z)

c11=z(1);
c12=z(2);
c13=z(3);
c33=z(4);
c44=z(5);
c66=z(6);
% Matriz de rigidez
rigidez=zeros(6,6);
rigidez(1,1)=c11;
rigidez(1,2)=c12;
rigidez(1,3)=c13;
rigidez(2,1)=c12;
rigidez(2,2)=c11;
rigidez(2,3)=c13;
rigidez(3,1)=c13;
rigidez(3,2)=c13;
rigidez(3,3)=c33;
rigidez(4,4)=c44;
rigidez(5,5)=c44;
rigidez(6,6)=c66;


compliance=inv(rigidez);

Et=1/compliance(1,1); % Et efectivo
Ea=1/compliance(3,3); % Ea efectivo
Gt=1/compliance(6,6); %Gt efectivo
Ga=1/compliance(4,4); %Ga efectivo
Pt=-compliance(1,2)/compliance(1,1); % Poisson t efectivo
Pa=-compliance(1,3)/compliance(3,3); % Poisson a efectivo
P31=-compliance(1,3)/compliance(1,1); % Poisson 31 efectivo

Et1=4*(c11+c12)/2*c44/(c11+c12)


%******************************************************
EE=[Ea Et Ga Gt Pa Pt];
