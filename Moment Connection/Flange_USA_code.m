clc
clear
%fprintf('Ry=1.5 for pipe and rolled box \nRy=1.2 for other rolled shape\nRu=1.15 for built up member\n\n')
%Ry=input('Ry= ');
o='Rafter.txt';
oo='Rafter';
OUTPUT = fopen(o,'w');
fprintf(OUTPUT,'EndPlate Connection for Beam to Column And Splice Beam\r\n');
d=134;
bf=40;
tf=2;
tw=1.2;
Ry=1.15;
fprintf(OUTPUT,'d=%2.2f cm\r\nbf=%2.2f cm\r\n',d,bf);
fprintf(OUTPUT,'tf=%2.2f cm\r\ntw=%2.2f cm\r\n',tf,tw);
%Mp=input('Mp= ');
Mp=0;

Fy=2400;
if (Mp==0)
    %Zx=input('Zx= ');
    Zx=3710;
    Mp=Zx*Fy;
end
Mpr=1.1*Ry*Mp;
%fprintf(OUTPUT,'Mu=Mpr=1.1*Ry*Mp=%2.2f t.m\r\n',Mpr/100000);
%Sh=input('Sh= ');
Sh=0;

if (Sh==0)
%     Mu=Mpr;
      Mu=660000;
else
    %Lh=input('Lh= ');
    %qu=input('qu= ');
    %L=input('L= ');
    
    Lh=1650;
    qu=1.41*2.4+0.2*7.52;
    L=1650;
    
     Mu=Mpr+(2*Mpr/Lh+qu*Lh/2)*Sh+qu*Sh.^2/2;
end
%Mu FOR DEMAND DESIGN
% Mu=345*100000;
Mu=345*100000;
fprintf(OUTPUT,'Mu=Maximum base on Beam to column connection and beam spice\n');
fprintf(OUTPUT,'Based on AISC341-10 CH.E Part E1-6b\n from ELM ID=143 \nMu(Ru=1)=%2.2f t.m\r\n',Mu/100000);
%Vu=input('Vu= ');
% Vu=13000;
 Vu=33000;
fprintf(OUTPUT,'Vu(Ru=1)=%2.2f ton\r\n',Vu/1000);
if (Vu==0)
    if (Sh==0)
    %Lh=input('Lh= ');
    %qu=input('qu= ');
    %L=input('L= ');
    Lh=1650;
    qu=1.41*2.4+0.2*7.52;
    L=1650;
    
    end        
Vu=2*Mpr/Lh+qu*L/2;
fprintf(OUTPUT,'Vu=2*Mpr/Lh+qu*L/2=%2.2f ton\r\n',Vu/1000);
end

%b=input('b= ');
b=40;
%H=input('H= ');
H=160;
fprintf(OUTPUT,'Dimension of EndPlate is %2.0f*%2.0f cm \r\n',b,H);

A=b*H;

%Nb=input('Number of Bolt= '); all bolts
Nb=24;

%Fub=input('Fu for Bolt= ');
Fub=10000;

%db=input('db= ');
db=3.6;
Ab=0.25*pi()*db.^2;
Tb=0.55*Fub*Ab;
Pp=Nb*Tb;
fpi=Pp/A;
fti=Mu*H/2/(b*H.^3/12);

fprintf(OUTPUT,'Nb=%2.0f\r\nFub=%2.0f kg/cm2\r\ndb=%2.1f cm\r\nAb=%2.2f cm2\r\n',Nb,Fub,db,Ab);
fprintf(OUTPUT,'Tb=0.55*Fub*Ab=%.0f kg \r\n',Tb);
%force because of pretension of Bolts
fprintf(OUTPUT,'Pp=Nb*Tb=%.0f kg \r\n',Pp);
fprintf(OUTPUT,'fpi=Pp/A=%2.0f kg/cm2 Stress between Plates \r\n',fpi);
fprintf(OUTPUT,'fti=Mu*H/2/(b*H.^3/12)=%2.0f kg/cm2 \r\n',fti);

if fti<fpi
    fprintf(OUTPUT,'fti/fpi=%2.2f  Number of Bolts are OK\r\n',(fti/fpi));
else
    fprintf(OUTPUT,'fti/fpi=%2.2f  Number of Bolts are NOT OK\r\n',(fti/fpi));
end

%%
%IBolt=input('I Bolt (Sigma Ad2)= ');

if mod(Nb,4)==0
    s=zeros(Nb/4,1);
    S=zeros(Nb/4,1);
else
    NumberColbolt=2;
    s=zeros((Nb+NumberColbolt)/4,1);
    S=zeros((Nb+NumberColbolt)/4,1);
end
s(1)=6.5;
s(2:5)=13;
s(6)=14;

S(1)=s(1);
for i=2:length(S)
    S(i)=s(i)+S(i-1);
end    
IBolt=0;
for i=1:length(S)
    IBolt=IBolt+4*.25*pi()*db.^2*(S(i)).^2;
end
fprintf(OUTPUT,'I-Bolts=%d cm4\r\n',IBolt);
%%


%CBolt=input('Top Bolt= ');
CBolt=S(end);

fbt=Mu*CBolt/IBolt;
fn=0.75*.75*Fub;
fprintf(OUTPUT,'fbt=Mu*CBolt/IBolt=%2.0f kg/cm2 Stress in Bolts\r\n',fbt);
fprintf(OUTPUT,'fn=0.75*.75*Fub=%2.0f kg/cm2 \r\n',fn);

if fbt<fn
    fprintf(OUTPUT,'Bolts are OK with Ratio %2.2f \r\n',fbt/fn);
else
    fprintf(OUTPUT,'Bolts are NOT OK with Ratio %2.2f \r\n',fbt/fn);
end

nb=1;
miu=0.5;
fi=1;
Du=1.13;
hf=1;
ns=1;
Tu=fbt*Ab;
Ksc=1-Tu/(Du*Tb*nb);
fiRni=fi*Ksc*miu*Du*hf*Tb*ns;
fiVn=Nb*fiRni;
fprintf(OUTPUT,'nb=%2.2f\r\nmiu=%2.2f\r\nfi=%2.2f \r\n',nb,miu,fi);
fprintf(OUTPUT,'Du=%2.2f\r\nhf=%2.2f\r\nns=%2.2f \r\n',Du,hf,ns);
fprintf(OUTPUT,'Tu=fbt*Ab=%2.2f kg\t\r\nKsc=1-Tu/(Du*Tb*nb)=%2.2f\r\n',Tu,Ksc);
fprintf(OUTPUT,'fiRni=fi*Ksc*miu*Du*hf*Tb*ns=%2.0f kg\r\nfiVn=Nb*fiRni=%2.0f kg\r\n',fiRni,fiVn);
fprintf(OUTPUT,'Vu/fiVn=%2.2f Shear Ratio is OK\r\n',Vu/fiVn);

%krishnamurthy
Aw=d*tw;
Af=bf*tf;
Tbu=Mu/(d-tf);
%FROM Bolt AXE TO UPPER FLANGE PAGE 480 AZHARI
pf=5.5;
aw=1;
b_prm=pf-.25*db-aw;
%Ca=1.14 FOR ST37& 10.9
%Ca=1.13 FOR ST37& 8.8
ca=1.14;
cb=sqrt(bf/b);
alfa_m=ca*cb*(Af/Aw).^(1/3)*(b_prm/db).^(1/4);
Me=alfa_m*Tbu*b_prm/4;
tp=sqrt(4*Me/(0.9*b*Fy));

fprintf(OUTPUT,'krishnamurthy method : \r\n');
fprintf(OUTPUT,'Aw=d*tw=%2.2f cm2\r\nAf=bf*tf=%2.2f cm2\r\n',Aw,Af);
fprintf(OUTPUT,'Tbu=Mu/(d-tf)=%2.0f kg\r\npf=%2.2f cm\r\n',Tbu,pf);
fprintf(OUTPUT,'b_prm=pf-0.25*db-aw=%2.2f cm\r\nCa=%2.2f\r\n',b_prm,ca);
fprintf(OUTPUT,'cb=sqrt(bf/b)=%2.2f \r\nalfa_m=ca*cb*(Af/Aw).^(1/3)*(b_prm/db).^(1/4)=%2.2f \t \r\n',cb,alfa_m);
fprintf(OUTPUT,'Me=alfa_m*Tbu*b_prm/4=%2.0f kg.cm \r\ntp=sqrt(4*Me/(0.9*b*Fy))=%2.2f cm \t \r\n',Me,tp);

fprintf(OUTPUT,'Thickness of plate must bigger than %2.2f cm \r\n',tp);
tp=2;
bst=13;
tst=1;
fprintf(OUTPUT,'USE tp=%2.2f for Plate and add stiffner with',tp);
fprintf(OUTPUT,' bst=%2.1f & tst=%2.1f\t \r\n',bst,tst);
ybar=(tp*b+bst*tst)/2/b;
Zst=ybar.^2/2*b+(tp-ybar).^2/2*b+(bst/2+tp-ybar)*bst*tst;
fi_Mnst=0.9*Zst*Fy;
fprintf(OUTPUT,'fi_Mnst=0.9*Zst*Fy=%2.0f kg.cm\r\n',fi_Mnst);
fprintf(OUTPUT,'Ratio of Mu/Mn for stiff & plate is %2.2f \r\n',Me/fi_Mnst);