clearvars
clf
ii=9;
 lam=[1e-10, .1:.1:.7];
% lam=[0, .1:.1:.5 ]; 
% disp(' ejes no principales=ejes global Harald')
% %%%***********************
% load  S_60_62_66_69
% S=S_62_66_69(:,ii);
% delta1=2*ZZ;
% Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
% w1=1;
% w2=Periodos2;
%%%*************************
% load S_Golavchan_Rect_square
% S1=S_Golavchan110(:,ii);
% delta1=2*ZZ(ii);
% Periodos2=Xg(ii).*(cos(theta(ii))+1i.*sin(theta(ii)));
% w1=1;
% w2=Periodos2;

%%%%********************
%****************************
%  load Sk_Tk_30_mas_5_675_hasta_90_ENoP
% S=Sk(:,ii);
 %%%%%%%%%%%%%%%%%%%%%

% load S10_45_75_90_delta1_2; 

%  S1=S10_45_75_90_delta1_2(:,ii);
 

S1=[0.000000000000000 + 0.000000000000000i
  0.000000000000000 + 0.000000000000000i
  0.000000000000000 + 0.000000000000000i
  2.073251677455220 + 1.196992414076610i
  0.000000000000000 + 0.000000000000000i
  2.345504840359280 - 2.345504840359280i
  0.000000000000000 + 0.000000000000000i
  1.228106433734542 + 2.127142740330432i
  0.000000000000000 + 0.000000000000000i
  3.486533339020110 - 0.934213792574630i
  0.000000000000000 + 0.000000000000000i
  0.000000000000008 - 0.196530919938053i
  0.000000000000000 + 0.000000000000000i
  3.852326730675952 + 1.032227836465435i
  0.000000000000000 + 0.000000000000000i
  0.956473396321423 - 1.656660518516660i
  0.000000000000000 + 0.000000000000000i
  1.959430896058406 + 1.959430896058379i];
delta1=3.301242952310458 - 0.182224507021762i;
theta=1.308996938995747;
Periodos2=(cos(theta)+1i*sin(theta));
 w1=1;
 w2=Periodos2;  
delta2=delta1.*w2-2*pi*1i; 
% cm=2-0.3*1i;
% cf=1-8*1i;
% delta1(ii)=(1.73289 - 0.07888i);
%  delta2=delta1.*w2/2-pi*1i/2; 
cm=1;
cf=120;
cmf=1;
x1=cmf/cm;
x2=cf/cm;
%lam=[.3 .5 .7];
 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
% Simetría cuadrada. c1313 y c2323 no perfect
% Caso Elastico y Transversalmente isotrópico Conatcto No Ideal.
%*************************************************************
% load S10_45_75_90
% Sk=S10_45_75_90;

% load Sarccos025
% Sk=Sarccos025;

% ************************************************************************

%  V = abs(w1(ii))*abs(w2(ii))*sin(theta(ii)); % *** Volumen de la celda periodica *******

V = abs(w1)*abs(w2)*sin(theta); % *** Volumen de la celda periodica *******

%%%%%*** expresiones de A y B en el problema I23 ********
% H1=(conj(delta1(ii))*conj(w2(ii))-conj(delta2(ii))*conj(w1(ii)))/(w1(ii)*conj(w2(ii))-w2(ii)*conj(w1(ii)));
% H2=((delta1(ii))*conj(w2(ii))-(delta2(ii))*conj(w1(ii)))/(w1(ii)*conj(w2(ii))-w2(ii)*conj(w1(ii)));

% 
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************

for j=1:length(lam)
  V3=lam(j);
  V2=1e-12;
  r1 = sqrt(V*(V2+V3)/pi);  % *** Radio de la fibra más la mesophase ************************************
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
 J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
 
 CENo0(j,:)=coeficientes_efectivos(x1,x2,V3,V2,X1,J); 
 CENo1(j,:) =ahm3f_Paralelogramo_Cortas_2023_No1(x1,x2,V3,V2,S1,X1,r1,J);
 CENo2(j,:) =ahm3f_Paralelogramo_Cortas_2023_No2(x1,x2,V3,V2,S1,X1,r1,J); 
 CENo3(j,:) =ahm3f_Paralelogramo_Cortas_2023_No3(x1,x2,V3,V2,S1,X1,r1,J);
 CENo4(j,:) =ahm3f_Paralelogramo_Cortas_2023_No4(x1,x2,V3,V2,S1,X1,r1,J);
 
end
CENoG(:,:) = ahm3f_elastico_p13_p23_75grados_Agosto;

matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
mesofase=[ 'cmf= ',num2str(cf,5)];
kappa1=[ '$\kappa_1$','=',num2str(cmf/cm,5),];
kappa2=[ '$\kappa_2$','=',num2str(cf/cm,5)];
angulo=['angle of the periodic cell= ' num2str(theta*180/pi,3),'$^o$'];
Modulo=[ '|w2|= ',num2str(abs(w2),5)];
disp(Modulo)
disp(angulo)
disp(matriz)
disp(fibra)
disp(mesofase)


Tabla(:,:)=[CENo0(:,1), CENo1(:,1),CENo2(:,1),CENo3(:,1),CENo4(:,1),CENoG(:,1)];

 disp('Vf   	c1313    c2323     c1323     c2313 ')    
 disp(' -------------------------------------------')
 disp(num2str([lam' Tabla(:,:)],6))
 
 disp(' -------------------------------------------')
 Tabla(:,:)=[CENo0(:,2), CENo1(:,2),CENo2(:,2),CENo3(:,2),CENo4(:,2),CENoG(:,2)];
 disp('Vf   	c1313    c2323     c1323     c2313 ')    
 disp(' -------------------------------------------')
 disp(num2str([lam' Tabla(:,:)],6))
 
 Tabla(:,:)=[CENo0(:,3), CENo1(:,3),CENo2(:,3),CENo3(:,3),CENo4(:,3),CENoG(:,3)];
 disp(' -------------------------------------------')
 disp('Vf   	c1313    c2323     c1323     c2313 ')    
 disp(' -------------------------------------------')
 disp(num2str([lam' Tabla(:,:)],6))
 

figure(1)
for jj=1:3
subplot(1,3,jj)
hold on

plot(lam,CENoG(:,jj) ,'ko-','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                        'MarkerSize',6)


plot(lam,CENo4(:,jj) ,'b<-','LineWidth',1,...
                       'MarkerEdgeColor','b',...
                       'MarkerFaceColor','b',...
                        'MarkerSize',4)

plot(lam,CENo3(:,jj) ,'gs-','LineWidth',1,...
                       'MarkerEdgeColor','g',...
                       'MarkerFaceColor','g',...
                        'MarkerSize',4) 
plot(lam,CENo2(:,jj) ,'rs-','LineWidth',1,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                        'MarkerSize',4) 
plot(lam,CENo1(:,jj) ,'ys-','LineWidth',1,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                        'MarkerSize',4) 
plot(lam,CENo0(:,jj) ,'mo-','LineWidth',1,...
                       'MarkerEdgeColor','m',...
                       'MarkerFaceColor','m',...
                        'MarkerSize',4)

hold off

ylabel('$\kappa$','fontsize',12,'interpreter','latex')
xlabel('$V_2$','fontsize',12,'interpreter','latex')
legend('$No=G$','$No=4$','$No=3$','$No=2$','$No=1$','$Van Fo Fi$','fontsize',12,'interpreter','latex')
grid on
box on
text(.1,CENo3(3,jj),kappa1,'fontsize',12,'fontsize',12,'interpreter','latex')
text(.1,CENo3(4,jj),kappa2,'fontsize',12, 'fontsize',12,'interpreter','latex')
text(.1,CENo3(5,jj),angulo,'fontsize',12, 'fontsize',12,'interpreter','latex')
end





% Funciones auxiliares
%%% orden del sistema No=1
function RR1 =ahm3f_Paralelogramo_Cortas_2023_No1(x1,x2,V3,V2,S1,X1,r1,J)

S4=S1(4);S6=S1(6);

X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);


N1=sqrt(3)*r1^4*[real(S4),-imag(S4);-imag(S4), -real(S4)];
N2=X3*N1;
W=10*r1^6*X3*[real(S6),-imag(S6);-imag(S6), -real(S6)];
AP2=eye(2)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR1 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end

%%% orden del sistema No=2
function RR2 =ahm3f_Paralelogramo_Cortas_2023_No2(x1,x2,V3,V2,S1,X1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);

X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);

N1=[sqrt(3)*r1^4*real(S4),-sqrt(3)*r1^4*imag(S4),sqrt(5)*r1^6*real(S6),-sqrt(5)*r1^6*imag(S6);
    -sqrt(3)*r1^4*imag(S4), -sqrt(3)*r1^4*real(S4),-sqrt(5)*r1^6*imag(S6), -sqrt(5)*r1^6*real(S6)];
N2=[X3*sqrt(3)*r1^4*real(S4),-X3*sqrt(3)*r1^4*imag(S4),  X5*sqrt(5)*r1^6*real(S6),-X5*sqrt(5)*r1^6*imag(S6);
    -X3*sqrt(3)*r1^4*imag(S4), -X3*sqrt(3)*r1^4*real(S4),-X5*sqrt(5)*r1^6*imag(S6), -X5*sqrt(5)*r1^6*real(S6)];

N2=N2';

W=[10*r1^6*X3*real(S6),-10*r1^6*X3*imag(S6),7*sqrt(15)*r1^8*X3*real(S8),-7*sqrt(15)*r1^8*X3*imag(S8);...
    -10*r1^6*X3*imag(S6),-10*r1^6*X3*real(S6),-7*sqrt(15)*r1^8*X3*imag(S8),-7*sqrt(15)*r1^8*X3*real(S8);...
    7*sqrt(15)*r1^8*X5*real(S8),-7*sqrt(15)*r1^8*X5*imag(S8),126*r1^10*X5*real(S10),-126*r1^10*X5*imag(S10);...
    -7*sqrt(15)*r1^8*X5*imag(S8),-7*sqrt(15)*r1^8*X5*real(S8),-126*r1^10*X5*imag(S10),-126*r1^10*X5*real(S10)];
AP2=eye(4)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR2 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end

%%% orden del sistema No=3
function RR3 =ahm3f_Paralelogramo_Cortas_2023_No3(x1,x2,V3,V2,S1,X1,r1,J)
              

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);

X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);

N1=[sqrt(3)*r1^4*real(S4),-sqrt(3)*r1^4*imag(S4),sqrt(5)*r1^6*real(S6),-sqrt(5)*r1^6*imag(S6),    sqrt(7)*r1^8*real(S8),-sqrt(7)*r1^8*imag(S8);
    -sqrt(3)*r1^4*imag(S4), -sqrt(3)*r1^4*real(S4),-sqrt(5)*r1^6*imag(S6), -sqrt(5)*r1^6*real(S6),-sqrt(7)*r1^8*imag(S8), -sqrt(7)*r1^8*real(S8)];
N2=[X3*sqrt(3)*r1^4*real(S4),-X3*sqrt(3)*r1^4*imag(S4),  X5*sqrt(5)*r1^6*real(S6),-X5*sqrt(5)*r1^6*imag(S6),X7*sqrt(7)*r1^8*real(S8),-X7*sqrt(7)*r1^8*imag(S8);
    -X3*sqrt(3)*r1^4*imag(S4), -X3*sqrt(3)*r1^4*real(S4),-X5*sqrt(5)*r1^6*imag(S6), -X5*sqrt(5)*r1^6*real(S6),-X7*sqrt(7)*r1^8*imag(S8), -X7*sqrt(7)*r1^8*real(S8)];
N2=N2';

W=[10*r1^6*X3*real(S6),-10*r1^6*X3*imag(S6),7*sqrt(15)*r1^8*X3*real(S8),-7*sqrt(15)*r1^8*X3*imag(S8),12*sqrt(21)*r1^10*X3*real(S10),-12*sqrt(21)*r1^10*X3*imag(S10);...
    -10*r1^6*X3*imag(S6),-10*r1^6*X3*real(S6),-7*sqrt(15)*r1^8*X3*imag(S8),-7*sqrt(15)*r1^8*X3*real(S8),-12*sqrt(21)*r1^10*X3*imag(S10),-12*sqrt(21)*r1^10*X3*real(S10);...
    7*sqrt(15)*r1^8*X5*real(S8),-7*sqrt(15)*r1^8*X5*imag(S8),126*r1^10*X5*real(S10),-126*r1^10*X5*imag(S10),66*sqrt(35)*r1^12*X5*real(S12),-66*sqrt(35)*r1^12*X5*imag(S12);...
    -7*sqrt(15)*r1^8*X5*imag(S8),-7*sqrt(15)*r1^8*X5*real(S8),-126*r1^10*X5*imag(S10),-126*r1^10*X5*real(S10),-66*sqrt(35)*r1^12*X5*imag(S12),-66*sqrt(35)*r1^12*X5*real(S12);...
    12*sqrt(21)*r1^10*X7*real(S10),-12*sqrt(21)*r1^10*X7*imag(S10),66*sqrt(35)*r1^12*X7*real(S12),-66*sqrt(35)*r1^12*X7*imag(S12),1716*r1^14*X7*real(S14),-1716*r1^14*X7*imag(S14);...
    -12*sqrt(21)*r1^10*X7*imag(S10),-12*sqrt(21)*r1^10*X7*real(S10),-66*sqrt(35)*r1^12*X7*imag(S12),-66*sqrt(35)*r1^12*X7*real(S12),-1716*r1^14*X7*imag(S14),-1716*r1^14*X7*real(S14)];

AP2=eye(6)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR3 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end

%%% orden del sistema No=3
function RR3 =ahm3f_Paralelogramo_Cortas_2023_No4(x1,x2,V3,V2,S1,X1,r1,J)
              

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);

X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);


         

N1=[sqrt(3)*r1^4*real(S4),-sqrt(3)*r1^4*imag(S4),sqrt(5)*r1^6*real(S6),-sqrt(5)*r1^6*imag(S6), sqrt(7)*r1^8*real(S8),-sqrt(7)*r1^8*imag(S8),3*r1^10*real(S10),-3*r1^10*imag(S10);
    -sqrt(3)*r1^4*imag(S4), -sqrt(3)*r1^4*real(S4),-sqrt(5)*r1^6*imag(S6), -sqrt(5)*r1^6*real(S6),-sqrt(7)*r1^8*imag(S8), -sqrt(7)*r1^8*real(S8),-3*r1^10*imag(S10), -3*r1^10*real(S10)];
N2=[X3*sqrt(3)*r1^4*real(S4),-X3*sqrt(3)*r1^4*imag(S4),  X5*sqrt(5)*r1^6*real(S6),-X5*sqrt(5)*r1^6*imag(S6),X7*sqrt(7)*r1^8*real(S8),-X7*sqrt(7)*r1^8*imag(S8),3*X9*r1^10*real(S10),-3*X9*r1^10*imag(S10);
    -X3*sqrt(3)*r1^4*imag(S4), -X3*sqrt(3)*r1^4*real(S4),-X5*sqrt(5)*r1^6*imag(S6), -X5*sqrt(5)*r1^6*real(S6),-X7*sqrt(7)*r1^8*imag(S8), -X7*sqrt(7)*r1^8*real(S8),-3*X9*r1^10*imag(S10), -3*X9*r1^10*real(S10)];
N2=N2';

W=[10*r1^6*X3*real(S6),-10*r1^6*X3*imag(S6),7*sqrt(15)*r1^8*X3*real(S8),-7*sqrt(15)*r1^8*X3*imag(S8),12*sqrt(21)*r1^10*X3*real(S10),-12*sqrt(21)*r1^10*X3*imag(S10),55*sqrt(27)*r1^12*X3*real(S12),-55*sqrt(27)*r1^12*X3*imag(S12);...
    -10*r1^6*X3*imag(S6),-10*r1^6*X3*real(S6),-7*sqrt(15)*r1^8*X3*imag(S8),-7*sqrt(15)*r1^8*X3*real(S8),-12*sqrt(21)*r1^10*X3*imag(S10),-12*sqrt(21)*r1^10*X3*real(S10),-55*sqrt(27)*r1^12*X3*imag(S12),-55*sqrt(27)*r1^12*X3*real(S12);...
    7*sqrt(15)*r1^8*X5*real(S8),-7*sqrt(15)*r1^8*X5*imag(S8),126*r1^10*X5*real(S10),-126*r1^10*X5*imag(S10),66*sqrt(35)*r1^12*X5*real(S12),-66*sqrt(35)*r1^12*X5*imag(S12),429*sqrt(5)*r1^14*X5*real(S14),-429*sqrt(5)*r1^14*X5*imag(S14);...
    -7*sqrt(15)*r1^8*X5*imag(S8),-7*sqrt(15)*r1^8*X5*real(S8),-126*r1^10*X5*imag(S10),-126*r1^10*X5*real(S10),-66*sqrt(35)*r1^12*X5*imag(S12),-66*sqrt(35)*r1^12*X5*real(S12),-429*sqrt(5)*r1^14*X5*imag(S14),-429*sqrt(5)*r1^14*X5*real(S14);...
    12*sqrt(21)*r1^10*X7*real(S10),-12*sqrt(21)*r1^10*X7*imag(S10),66*sqrt(35)*r1^12*X7*real(S12),-66*sqrt(35)*r1^12*X7*imag(S12),1716*r1^14*X7*real(S14),-1716*r1^14*X7*imag(S14),2145*sqrt(7)*r1^16*X7*real(S16),-2145*sqrt(7)*r1^16*X7*imag(S16);...
    -12*sqrt(21)*r1^10*X7*imag(S10),-12*sqrt(21)*r1^10*X7*real(S10),-66*sqrt(35)*r1^12*X7*imag(S12),-66*sqrt(35)*r1^12*X7*real(S12),-1716*r1^14*X7*imag(S14),-1716*r1^14*X7*real(S14),-2145*sqrt(7)*r1^16*X7*imag(S16),-2145*sqrt(7)*r1^16*X7*real(S16);...
      55*sqrt(27)*r1^12*X9*real(S12), -55*sqrt(27)*r1^12*X9*imag(S12),429*sqrt(5)*r1^14*X9*real(S14),-429*sqrt(5)*r1^14*X9*imag(S14), 2145*sqrt(7)*r1^16*X9*real(S16),-2145*sqrt(7)*r1^16*X9*imag(S16),24310*r1^18*X9*real(S18),-24310*r1^18*X9*imag(S18);...
      -55*sqrt(27)*r1^12*X9*imag(S12), -55*sqrt(27)*r1^12*X9*real(S12),-429*sqrt(5)*r1^14*X9*imag(S14),-429*sqrt(5)*r1^14*X9*real(S14),-2145*sqrt(7)*r1^16*X9*imag(S16),-2145*sqrt(7)*r1^16*X9*real(S16),-24310*r1^18*X9*imag(S18),-24310*r1^18*X9*real(S18)];
     
AP2=eye(8)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR3 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end

function RRG=ahm3f_elastico_p13_p23_75grados_Agosto
RRG=[1.000000000196694   1.000000000196694   0.000000000000000   0.000000000000000
   1.218518065366796   1.217803588398004   0.001333232174188   0.001333232174188
   1.491758861903358   1.488156265118038   0.006722537120996   0.006722537120996
   1.844284177641961   1.833740316476588   0.019675112788564   0.019675112788564
   2.320161785953333   2.294765272644918   0.047390539001056   0.047390539001056
   3.008842268195630   2.952026083620859   0.106020443762628   0.106020443762628
   4.128840818021787   4.001232525925717   0.238120314784820   0.238120314784820
   6.437477626338564   6.115340106509423   0.601116795513297   0.601116795513297];

end

%**************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
%**************************************************************************
function RR =coeficientes_efectivos(x1,x2,V3,V2,X1,Z)
F=X1*(1-x1+2*V3*(x1-x2)*x1/(V2*(x1+x2)+2*V3*x1));
G=1/(V2*(x1+x2)+2*V3*x1);
DtZ=det(Z);
c55=(1-V2-V3+x1*V2+x2*V3-(V2+V3)*F*((X1+1)*Z(2,2)-DtZ)/(DtZ*X1)-G*(x1-x2)^2*V2*V3);
                    c45=((V2+V3)*F*((X1+1)*Z(2,1))/(DtZ*X1));
                    c54=((V2+V3)*F*((X1+1)*Z(1,2))/(DtZ*X1));
c44=(1-V2-V3+x1*V2+x2*V3-(V2+V3)*F*((X1+1)*Z(1,1)-DtZ)/(DtZ*X1)-G*(x1-x2)^2*V2*V3);
      
RR=[c55,c44,c54,c45];
end
