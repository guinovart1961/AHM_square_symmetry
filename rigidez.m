function ce=rigidez(Z)

Ea=Z(1);
Et=Z(2);
Ga=Z(3);
Gt=Z(4);
Pa=Z(5);
Pt=Z(6);

mat=zeros(6,6);
mat(1,1)=1/Et;
mat(1,2)=-Pt/Et;
mat(1,3)=-Pa/Ea;
mat(2,1)=-Pt/Et;
mat(2,2)=1/Et;
mat(2,3)=-Pa/Ea;
mat(3,1)=-Pa/Ea;
mat(3,2)=-Pa/Ea;
mat(3,3)=1/Ea;
mat(4,4)=1/Ga;
mat(5,5)=1/Ga;
mat(6,6)=1/Gt;


nat=inv(mat);


c11=nat(1,1);
c12=nat(1,2);
c13=nat(1,3);
c33=nat(3,3);
c44=nat(4,4);
c66=nat(6,6);
k=(c11+c12)/2;
ce=[c11 c12 c13 c33 c44 c66 k];