%Resuelve el problema elástico bifásico fibroso con simetria cuadrada
%      

clear
ang=9; % 90
load Sk_Tk_304045506070758090_ENoP_E

% Tsai y Hann 1980
%epoxy                          Glass
Ea1= 3.45;                 Ea2= 73.1;                   
Et1= 3.45;                 Et2= 73.1;                    
Ga1= 1.2778;               Gt2= 29.9590;               
Pa1= 0.35 ;                Pa2= 0.22;                                    
Pt1= 0.35;                 Pt2= 0.22;                    
Gt1=Et1/(2*(1+Pt1)); Ga2=Ea1/(2*(1+Pa1)); 

%*************************************************************
%INTRODUCCION DE LOS DATOS DE LA MATRIZ                      *
%*************************************************************
Z=[Ea1,Et1,Ga1,Gt1,Pa1,Pt1];
Im=Z;
cm=rigidez(Z); % coeficientes elásticos 
%*************************************************************
%INTRODUCCION DE LOS DATOS DE LA FIBRA                       *
%*************************************************************

Z=[Ea2,Et2,Ga2,Gt2,Pa2,Pt2];
If=Z;
cf=rigidez(Z); % coeficientes elásticos 
%************************************************************
lam2=[0:.01:.1, .15:.05:.75, .78 .784]; %%%input('Fracción Volumétrica de la Fibra = ');
% disp('Tipo de simetria de las fibras: cuadrada(=1), hexagonal(<>1) ')
% sim=input('Tipo se simetria = ');

%************************************************************
% Ciclo para las diferentes combinaciones de fibra y mesofase
for j=1:length(lam2)
%**********************************************
%RR=AHMelasticoSq(cm,cf,lam2(j));
RR=AHMelasticoSqr_2010(cm,cf,lam2(j),Sk(:,ang),Tk(:,ang),gam1(ang),delta1(ang),theta(ang));
CE(j,:)=[RR];
CE1(j,:)=[lam2(j) RR]; % para la salida de los coeficietes

end


%*******************************************************
%*******************************************************
% Para hallar las Constantes Ingenieriles Efectivas

for j=1:length(lam2)
% EE1=ingenieriles3_Hashin79(CE(j,:)); % formulas en funcion de las contantes elasticas 
EE1=ingenieriles3(CE(j,:)); % 
EE(j,:)=EE1; %Constantes Ingenieriles Efectivas
end

%***************************************

% Para hallar las Cotas de Hashin 79

for j=1:length(lam2)
[LB,UB]=Hashin79(Im,If,cm,cf,lam2(j));
CHE(j,:)=[lam2(j) LB(1) EE(j,1) UB(1)  LB(2) EE(j,2) UB(2)]; %Cotas de Hashin (Modulos de Young)
CHG(j,:)=[lam2(j) LB(3) EE(j,3) UB(3)  LB(4) EE(j,4) UB(4)]; %Cotas de Hashin (Shears)
CHP(j,:)=[lam2(j) LB(5) EE(j,5) UB(5)  LB(6) EE(j,6) UB(6)]; %Cotas de Hashin (Poison)

end

%***************************************
% Salida
disp('         ***********************************************')
disp('            Coeficientes Efectivos Square Symmetry')
disp('***********************************************************************')
disp('    Vf        C11       C12         C13        C33       C44         C66 ')
disp('***********************************************************************')
disp(CE1)


% Salida
disp('         ***********************************************')
disp('            Cotas de Hashin 1979 (Modulos de Young)')
disp('***********************************************************************')
disp('    Vf        LEa       Ea         UEa        LEt       Et         UEt ')
disp('***********************************************************************')
disp(CHE)


% Salida
disp('         ***********************************************')
disp('            Cotas de Hashin 1979 (Modulos de Cizalladura)')
disp('***********************************************************************')
disp('    Vf        LGa       Ga         UGa        LGt       Gt         UGt ')
disp('***********************************************************************')
disp(CHG)

% Salida
disp('         ***********************************************')
disp('            Cotas de Hashin 1979 (Modulos de Poison)')
disp('***********************************************************************')
disp('    Vf        LPa       Pa         UPa        LPt       Pt         UPt ')
disp('***********************************************************************')
disp(CHP)

 h=figure;
  %*********** Comparación E_{1} 
  subplot(3,2,1)
plot(lam2,CHE(:,4)/Im(1),'k-.',lam2,CHE(:,3)/Im(1),'r-',...
    lam2,CHE(:,2)/Im(1),'k:')
ylabel('E_{a}/E_1(GPa)')
% axis([0 .8  0  380])
legend('CSH','MHA','CIH')

 subplot(3,2,2)
plot(lam2,CHE(:,7)/Im(2),'k-.',lam2,CHE(:,6)/Im(2),'r-',...
    lam2,CHE(:,5)/Im(2),'k:')
ylabel('E_{t}/E_1(GPa)')
%axis([0 .7  0  150])


  %*********** Comparación G 
  subplot(3,2,3)
plot(lam2,CHG(:,4)/Im(3),'k-.',lam2,CHG(:,3)/Im(3),'r-',...
    lam2,CHG(:,2)/Im(3),'k:')
ylabel('G_{a}/G_1(GPa)')
%axis([0 .7  0  150])

 subplot(3,2,4)
plot(lam2,CHG(:,7)/Im(4),'k-.',lam2,CHG(:,6)/Im(4),'r-',lam2,CE(:,7)/Im(4),'r:',...
    lam2,CHG(:,5)/Im(4),'k:')
ylabel('G_{t}/G_1, m^{*}/G_1(GPa)')
%axis([0 .7  0  150])

  %*********** Comparación \nu 
  subplot(3,2,5)
plot(lam2,CHP(:,4)/Im(5),'k-.',lam2,CHP(:,3)/Im(5),'r-',...
    lam2,CHP(:,2)/Im(5),'k:')
xlabel('Volumen de fibra (V_{2})') ,ylabel('\nu_{a}/\nu_{1}')
%axis([0 .7  0  150])

 subplot(3,2,6)
plot(lam2,CHP(:,7)/Im(6),'k-.',lam2,CHP(:,6)/Im(6),'r-',...
    lam2,CHP(:,5)/Im(6),'k:')
xlabel('Volumen de fibra (V_{2})') ,ylabel('\nu_{t}/\nu_{1}')
%axis([0 .7  0  150])
