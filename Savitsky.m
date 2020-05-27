% implementing http://www.westlawn.edu/ReferenceInfo/SavitskyPlaningHulls1964.pdf
function Savitsky
  clc

  ##g=9.81
  ##rho=1025
  ##nu=1.787e-6;
  ##beta=10 %grader
  ##V=40*1852/3600
  ##b=14*0.3048
  ##depl=60000*0.453592*g % i N
  ##LCG=29*0.3048
  ##VCG=2*0.3048
  ##epsilon = 4
  ##f=0.5*0.3048

  g=32.17405
  rho=1.95
  nu=1e-5
  beta=10 %grader
  V=67.5
  b=14
  depl=60000
  LCG=29
  VCG=2
  epsilon = 4
  f=0.5

  tau=3
  tau11=power(tau,1.1)


  CV=V/sqrt(g*b)

  CLbeta=depl/(0.5*rho*power(V,2)*power(b,2))

  %with deadrise

  CL0=findx(@CLbetaCalc, CLbeta, beta)

  CL0divtau11=CL0/tau11

  % zero deadrise:
  %CL0=tau11*(0.012*power(lambda, 0.5) + 0.0055*power(lambda, 2.5)/power(Cv,2)) % (15)


  % 0.60 <= Cv <= 13.00, 2 <= tau <= 15, lambda <= 4

  %%%lambda=(d/sin(tau) - b*tand(beta)/(2*pi*tand(tau)))/b % wetted length-beam ratio
  lambda=findx(@CL0divtau11Calc, CL0divtau11, CV)

  V1=V*sqrt(1-0.012*power(tau, 1.1)/(power(lambda, 0.5)*cosd(tau)))


  Re=V1*lambda*b/nu

  Cf=0.00174
  Cf=0.00184
  deltaCf=0.0004

  ##lambda=3.85

  Df=rho*V1**2*lambda*b**2*(Cf+deltaCf)/(2*cosd(beta))

  ##DfLbs=(Df/9.81)/0.453592

  tand(tau)
  %row 12
  sind(tau)
  cosd(tau)

  % row 14:
  row14=depl*tand(tau)

  % row 15:
  row15=Df/cosd(tau)
  %row 16:
  D=depl*tand(tau)+Df/cosd(tau)

  % row 17:
  Cp=0.75 - 1/((5.21*CV**2)/lambda**2+2.39)

  %row 18
  row18=Cp*lambda*b
  %row 19
  c=LCG-Cp*lambda*b

  %row 21
  a=VCG-(b/4)*tand(beta)

  %row 22
  row22=sind(tau+epsilon)
  %row 23
  row23=1-sind(tau)*sind(tau+epsilon)

  %row 24
  row24=row23*(c/cosd(tau))
  row24=(1-sind(tau)*sind(tau+epsilon))*(c/cosd(tau))

  row25=f*sind(tau)

  row26=row24-row25

  row27=depl*row26

  row28=a-f
  row29=Df*(a-f)

  row30=row27+row29

  N=depl*(1-sind(tau)*sind(tau+epsilon))/cosd(tau) % (34)

  eq35=depl*( (1-sind(tau)*sind(tau+epsilon))*c/cosd(tau) - f*sind(tau)) + Df*(a-f) % should be =0


end

function x=findx(func, yGoal, args)
  x0=1;
  y0=func(x0, args);
  
  x=x0;

  dx=1/10000000;
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


function CLbeta=CLbetaCalc(CL0, args)
  beta=args(1);
  CLbeta=CL0-0.0065*beta*power(CL0, 0.6); % (16) figure 11 
end

function CL0divtau11=CL0divtau11Calc(lambda, args)
  CV=args(1);
  CL0divtau11=(0.012*power(lambda, 0.5) ...
    + 0.0055*power(lambda, 2.5)/power(CV,2));
end







