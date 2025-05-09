%% Matlab file 
clc
clear


OUTPUT = fopen('BasePlate.txt','w');
fprintf(OUTPUT,'BasePlate Connection\r\n');
%fprintf('Ry=1.5 for pipe and rolled box \nRy=1.2 for other rolled shape\nRu=1.15 for built up member\n\n')
Ry=1.15;
d=37.4;
bf=40;
tf=2;
tw=1.2;
Fy=2400;
n=10;
fc=250;
B=80;
D=80;
A1=B*D;
A2=100*100;
f=30;
As=12*3.14*3*3/4;
Fub=6000;
Z=3327.5;
H=355;

fprintf(OUTPUT,'d=%2.2f cm\r\nbf=%2.2f cm\r\n',d,bf);
fprintf(OUTPUT,'tf=%2.2f cm\r\ntw=%2.2f cm\r\n',tf,tw);
fprintf(OUTPUT,'Ry=%2.2f\r\nFy=%2.0f kg/cm2\r\n',Ry,Fy);
fprintf(OUTPUT,'n=%2.0f\tfc=%2.0f kg/cm2\t',n,fc);
fprintf(OUTPUT,'B=%2.0f cm\tD=%2.0f cm\n',B,D);
fprintf(OUTPUT,'A1=%2.0f cm2\tA2=%2.0f cm2\t',A1,A2);
fprintf(OUTPUT,'f=%2.2f cm\n',f);
fprintf(OUTPUT,'As=%2.2f cm2\tFub=%2.0f cm\n',As,Fub);


for i=1:1
    if i==1
%MAX FORCE FORM COMBO
fprintf(OUTPUT,'Pu :Based on Amplified seismic load Combinations \n');
fprintf(OUTPUT,'Mu :Based on min(Amplified seismic load ,1.2RyFyZ )\n');
fprintf(OUTPUT,'Vu :Based on min(Amplified seismic load ,2RyFyZ/H )\n');

Mu=min(42*100000,1.1*Ry*Z*Fy);
Vu=min(16*1000,2*Ry*Fy*Z/H);
Pu=96*1000;
fprintf(OUTPUT,'Mu=%2.2f t.m\t',Mu/100000);
fprintf(OUTPUT,'Vu=%2.2f ton\t',Vu/1000);
fprintf(OUTPUT,'Pu=%2.2f ton\n',Pu/1000);

    elseif i==2
%MAX FORCE FORM COMBO OF OMEGA0
fprintf(OUTPUT,'\nResult for Amplified seismic load Combinations\n');
Mu=45*100000;
Vu=15*1000;
Pu=90*1000;
fprintf(OUTPUT,'Pu=%2.0f t\r\n',Pu/1000);
    else
%MAX FORCE MP
fprintf(OUTPUT,'\nResult for Plastic moment\n');
Mu=1.1*Ry*3327*Fy;
Vu=10;
Pu=10;
fprintf(OUTPUT,'Mu=Mpr=1.1*Ry*Mp=%2.0f t.m\r\n',Mu/100000);
    end


%Mp=input('Mp= ');
% Mp=0;

% if (Mp==0)
%     %Zx=input('Zx= ');
%     Zx=285;
%     
%     Mp=Zx*Fy;
% end
% Mpr=1.1*Ry*Mp;
% %n=Es/Ec

e=Mu/Pu;
fprintf(OUTPUT,'e=Mu/Pu=%2.0f cm\n',e);

a1=3*(e-B/2);
a2=6*n*As/D*(f+e);
a3=-a2*(B*.5+f);

% Coefficients of the cubic
coeffs = [1, a1, a2, a3];
r = roots(coeffs);            % find all (possibly complex) roots
real_roots = r(abs(imag(r))<1e-6);  % keep only nearly-real ones
if isempty(real_roots)
    error('No real root found');
end
x = max(real_roots);          % or choose the physically meaningful root

   
fup=2*Pu*(e+f)/(x*D*(B/2+f-x/3));
Tu=Pu*(e+x/3-B/2)/(B/2+f-x/3);
fiFp=0.65*min(0.85*fc*sqrt(A2/A1),1.7*fc);

fprintf(OUTPUT,'x=%2.2f \n',x);
fprintf(OUTPUT,'fup=2*Pu*(e+f)/(x*D*(B/2+f-x/3))=%2.2f\n',fup);
fprintf(OUTPUT,'Tu=Pu*(e+x/3-B/2)/(B/2+f-x/3)=%2.2f\n',Tu);
fprintf(OUTPUT,'fiFp=0.65*min(0.85*fc*sqrt(A2/A1),1.7*fc)=%2.2f\n',fiFp);

if fup<fiFp
    fprintf(OUTPUT,'fup/fiFp=%2.2f so its OK\n',fup/fiFp);
else
    fprintf(OUTPUT,'****fup/fiFp is NG (%2.2f)***** \n',fup/fiFp);
end
fut=Tu/As;
fuv=Vu/2/As;
Fnt=0.75*Fub;Fnv=0.45*Fub;
Fpnt=min(Fnt*(1.3-fuv/0.75/Fnv),Fnt);
Fpnv=min(Fnv*(1.3-fut/0.75/Fnt),Fnv);

fprintf(OUTPUT,'fut=Tu/As=%2.2f\n',fut);
fprintf(OUTPUT,'fuv=Vu/2/As=%2.2f\n',fuv);
fprintf(OUTPUT,'Fnt=0.75*Fub=%2.2f\n',Fnt);
fprintf(OUTPUT,'Fnv=0.45*Fub=%2.2f\n',Fnv);
fprintf(OUTPUT,'Fpnt=min(Fnt*(1.3-fuv/0.75/Fnv),Fnt)=%2.2f\n',Fpnt);
fprintf(OUTPUT,'Fpnv=min(Fnv*(1.3-fut/0.75/Fnt),Fnv)=%2.2f\n',Fpnv);


if fut<0.75*Fpnt
    fprintf(OUTPUT,'fut/0.75*Fpnt=%2.2f so its OK\n',fut/(0.75*Fpnt));
else
    fprintf(OUTPUT,'****fut/0.75*Fpnt is NG (%2.2f)**** \n',fut/(0.75*Fpnt));
end
if fuv<0.75*Fpnv
    fprintf(OUTPUT,'fuv/0.75*Fpnv=%2.2f so its OK\n',fuv/(0.75*Fpnv));
else
    fprintf(OUTPUT,'****fut/0.75*Fpnv is NG (%2.2f)**** \n',fuv/(0.75*Fpnv));
end
%for I
m=(B-.95*d)/2;
n=(D-0.8*bf)/2;
Mu1p=fup*m^2/3+fup*(x-m)/x*m.^2/6;
Mu2p=(fup+fup*(x-m)/x)/2*n.^2/2;
Mup=max(Mu1p,Mu2p);
tp_nostiff=2.11*sqrt(Mup/Fy);

fprintf(OUTPUT,'m=(B-.95*d)/2=%2.2f\n',m);
fprintf(OUTPUT,'n=(D-0.8*bf)/2=%2.2f\n',n);
fprintf(OUTPUT,'Mu1p=fup*m^2/3+fup*(x-m)/x*m.^2/6=%2.2f\n',Mu1p);
fprintf(OUTPUT,'Mu2p=(fup+fup*(x-m)/x)/2*n.^2/2=%2.2f\n',Mu2p);
fprintf(OUTPUT,'Mup=max(Mu1p,Mu2p)=%2.2f\n',Mup);
fprintf(OUTPUT,'thickness of plate without stiffener is %2.2f\n',tp_nostiff);

b=B;
tp=2.5;
bst=20.5;
tst=1.2;
numStiff=3;
fprintf(OUTPUT,'Try stiffener:tp=%2.2f\tbst=%2.2f\ttst=%2.2f\n',tp,bst,tst);

%Three STIFF
% ybar=(tp*b+3*bst*tst)/2/b;
% Zst=ybar.^2/2*b+(tp-ybar).^2/2*b+(bst/2+tp-ybar)*bst*tst*3;

%TWO STIFF
% ybar=(tp*b+2*bst*tst)/2/b;
% Zst=ybar.^2/2*b+(tp-ybar).^2/2*b+(bst/2+tp-ybar)*bst*tst*2;

%n STIFF
ybar=(tp*b+bst*tst*numStiff)/2/b;
Zst=(ybar.^2/2*b+(tp-ybar).^2/2*b+(bst/2+tp-ybar)*bst*tst*numStiff);

fiMnst=0.9*Zst*Fy;
Must=b*Mup;
% fprintf('thikness of stiffner is %2.2f \n',tst)
fprintf(OUTPUT,'Mu/Mn for stiffener & plate is %2.2f So its OK\n',Must/fiMnst);
end