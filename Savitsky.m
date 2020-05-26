function Savitsky
clc
%clear all
%# tau trim
%
%b=2.3 # beam
%beta=15 # deamrise
%tau=5 # angle of attack
%d=0.4 # draft at speed
%
%L2=b *tand(beta)/(2*tand(tau)) # water level line length
%
%lambda=(d/sin(tau) - b*tand(beta)/(2*pi*tand(tau)))/b % wetted length-beam ratio
%
%Lk=d/sin(tau)
%lambda=(Lk+Lc)/(2*b)
%
%# phi is the angle between the keel and spray edge measured in the plane of the bottom
%# tand(phi)=A+k1/(1-Ak1)
%A=((sind(tau)²*(1-2*K)+K²tand(tau)²*power(1/sind(beta)²-sind(tau)²),1/2))/(cosd(tau)+K*tand(tau)*sind(tau))
%k1=K*tand(tau)/sind(beta)
%
%tand(theta)=tand(phi)*cosd(beta)
%
%tand(alpha)=(pi/2)*(tand(tau)/tand(beta))
%
%As=(b²/2)*(tand(beta)/(pi*tand(tau)) - 1/(4*tand(phi)*cosd(beta)))
%
%CL=A*tau+B*tau²
%
%CLd=c*power(lambda, 0.5)*power(tau, 1.1) % 10

% Hydrostatic lift:
%Lb=(1/2)*rho*g*b³*(lambda-0.3)²*tand(tau) # 11
%
%CL=power(tau,1.1)*(c*power(lambda, 0.5) + D*power(lambda, n)/Cv²) % 14
g=9.81
rho=1025
beta=10 %grader
V=40*1852/3600
b=14*0.3048
depl=60000*0.453592*g % i N
LCG=29*0.3048
VCG=2*0.3048
epsilon = 4
f=0.5*0.3048

CV=V/sqrt(g*b)

CLbeta=depl/(0.5*rho*power(V,2)*power(b,2))

%with deadrise

CL0=findx(@CLbetaCalc, CLbeta, beta)

%CLbetaCalc=0;
%CL0=0;
%while(abs(CLbetaCalc-CLbeta)>0.00001)
%  CL0=CL0+0.000001;
%  CLbetaCalc=CL0-0.0065*beta*power(CL0, 0.6); % (16) figure 11 
%end
%CL0

tau=2
tau11=power(tau,1.1)

CL0divtau11=CL0/tau11

% zero deadrise:
%CL0=tau11*(0.012*power(lambda, 0.5) + 0.0055*power(lambda, 2.5)/power(Cv,2)) % (15)


% 0.60 <= Cv <= 13.00, 2 <= tau <= 15, lambda <= 4

%%%lambda=(d/sin(tau) - b*tand(beta)/(2*pi*tand(tau)))/b % wetted length-beam ratio
lambda=findx(@CL0divtau11Calc, CL0divtau11, CV)

V1=V*sqrt(1-0.012*power(tau, 1.1)/(power(lambda, 0.5)*cosd(tau)))

nu=1.787e-6;

Re=V1*lambda*b/nu

Cf=0.00174;
deltaCf=0.0004;

Df=rho*V**2*lambda*b**2*(Cf+deltaCf)/(2*cosd(beta))

DfLbs=(Df/9.81)/0.453592

tand(tau)
%row 12
sind(tau)
cosd(tau)

% row 14:
depl*tand(tau)

% row 15:
Df/cosd(tau);
%row 16:
D=depl*tand(tau)+Df/cosd(tau)

% row 17:
Cp=0.75 - 1/((5.21*CV**2)/lambda**2+2.39)

%row 18
Cp*lambda*b;
%row 19
c=LCG-Cp*lambda*b

%row 21
a=VCG-(b/4)*tand(beta)

%row 22
row22=sind(tau+epsilon)
%row 23
row23=1-sind(tau)*sind(tau+epsilon)

%row 24
row24=(1-sind(tau)*sind(tau+epsilon))*(c/cosd(tau))

row25=f*sind(tau)

row26=row24+row25

row27=depl*row26

row28=a-f
row29=Df*(a-f)

row30=row27+row29

eq35=depl*(row

% zero deadrise:
%CL0=tau11*(0.012*power(lambda, 0.5) + 0.0055*power(lambda, 2.5)/power(Cv,2)) % 15
%CL0divtau11=(0.012*power(lambda, 0.5) + 0.0055*power(lambda, 2.5)/power(Cv,2))
% 0.60 <= Cv <= 13.00, 2 <= tau <= 15, lambda <= 4
%
% Dp=delta*tand(tau) % 17
%
% D=delta*tand(tau) + Df/cosd(tau) % 18
%
% Df=Cfpho*V1²*(lambda*b²)/(2*cosd(beta)⁴) % 19
%
% CLd=0.012*power(lambda, 0.5)*power(tau, 1.1) % 20
%
% deltad=(1/2)*rho*V²*b²*(0.012*power(lambda, 0.5)*power(tau, 1.1)) % 20 dynamic loadon the bottom
%
% pd=deltad/(lambda*b²*cosd(tau))
%
% V1=V*power(1-2*pd/(rho*V²), 0.5) % (23) average bottom velocity
%
% trim larger than 4 deg (9) and (17):
% D=delta*tand(tau) + rho*V1²*Cf*lambda*b²/(2*cosd(beta)*cosd(tau))
%
% Rn=V1*lambda*b/nu # where nu is the kinematic viscosity
%
% D/delta=tand(tau)+(rho*V1²*Cf*lambda*b²)/(2*delta*cosd(beta)*cosd(tau)) # (26) drag-lift ratio
%
% D/delta=tand(tau) + ((V1/V)²*Cf*lambda)/(CL*cosd(tau)*cosd(beta)) % (27)

% vertical equilibrium of forces:
% delta0=N*cosd(tau) + T*sind(tau+epsilon)-Df*sind(tau) % (29)

# for horizontal equilibrium of forces:

%%T*cosd(tau+epsilon)=Df*cosd(tau) + N*sind(tau) % (30)

% for equilibrium of pitching moments:
%%Nc + Df*a - Tf=0 % (31)

% equilibrium trim taue
% trim at which (30) = 0



%%%CV=V/sqrt(g*b)
%CLbeta= delta/((1/2)*rho*V²*b²)


end

function x=findx(func, yGoal, args)
  
  x0=1;
  y0=func(x0, args);
  
  x=x0;

  dx=1/100000;
  y=y0;
  y=func(x, args);
  
  while(abs(yGoal-y)>0.00000001)
    dy=func(x+dx, args)-func(x, args);
    deriv=dx/dy;
    
    yGoal;
    x=x+(yGoal-y)*deriv;
    y=func(x, args);
    
  end
  
end

function y=funfun(x, args)
  y=x**2;
end

%lambda=0;
%CL0divtau11Calc=0;
%while(abs(CL0divtau11Calc-CL0divtau11)>0.0001)
%    
%  dlambda=
%  
%  lambda=lambda+0.000001;
%  CL0divtau11Calc=(0.012*power(lambda, 0.5) + 0.0055*power(lambda, 2.5)/power(CV,2));
%end

function CLbeta=CLbetaCalc(CL0, args)
  beta=args(1);
  CLbeta=CL0-0.0065*beta*power(CL0, 0.6); % (16) figure 11 
end

function CL0divtau11=CL0divtau11Calc(lambda, args)
  CV=args(1);
  CL0divtau11=(0.012*power(lambda, 0.5) ...
    + 0.0055*power(lambda, 2.5)/power(CV,2));
end







