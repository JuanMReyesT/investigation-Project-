clear all
WT.nB=3;
WT.R=0.61+0.08;
xs=0.160;
xt=WT.R;
rstep=0.04;
Total_Power=0;
Total_Thrust=0;
Total_Torque=0;
time_steps=120;
t=0;
%Wind turbine characteristics 
WT.rs_in1=[0;0;0.17018];
WT.rt_in1=[4;0;0];
WT.r=[xs:rstep:xt];
WT.chord=[0.1423 0.127 0.106 0.0875 0.0715 0.0575 0.0466 0.0402 0.037 0.035 0.0327 0.0304 0.028 0.02571]; %chord (m)
WT.twist=[27.03 28.68 27.03 23.93 22.13 19.73 18.33 18.33 18 16 14.93 13.53 12.73 12.03 11.5]*pi/180;
WT.dr=rstep ;
WT.rhub=0.08;
WT.ty=0;
WT.tt=0;
WT.tw=0;
WT.tc=0;


Sim.rho=1.225 ;
Sim.KinVisc=1.470e-5 ; % Kinematic viscosity [m^2/s] (for Reynolds number)
Sim.PITCH=0.0 ; % Pitch angle [rad] |+
Sim.dt=0.01;%delta time
Sim.RPM=650 ; % Rotational velocity [rad/s]
Wind.V0=12.5;% wind velocity
At=Sim.dt;
psi=0;
Omega=Sim.RPM*2* pi /60;
x=[0,psi];
v=[0,Sim.RPM*2* pi /60 ];
%algorithm parameters 
Algo.nbIt=60 ; % Maximum number of iterations
Algo.aTol= 0.01; % Tolerance in axial induction
Algo.bTipLoss=true; % True if tip -losses are applied
%----initialization of NMO
NMO.W_qs=ones(3,length(WT.r),WT.nB);
NMO.W0=ones(3,length(WT.r),WT.nB);
NMO.W=ones(3,length(WT.r),WT.nB);
NMO.W_int=ones(3,length(WT.r),WT.nB);
NMO.khi=0;
NMO.fs=ones(3,length(WT.r),WT.nB);
Power_vectorY=zeros(1,time_steps);
Power_thrustY=zeros(1,time_steps);
Power_TorqueY=zeros(1,time_steps);
Time=zeros(1,time_steps);
% % init W factors
% % initialization 
for o=1:42
[RES,NMO]=fBEMUnsteadyYaw(WT,Sim,Wind,Algo,NMO,x,v);
psi=psi+Omega*At;
x=[0,psi];
end
psi=0;
x=[0,psi];
%Simulation in time
for o=1:time_steps
Time(o)=t;
[RES,NMO]=fBEMUnsteadyYaw(WT,Sim,Wind,Algo,NMO,x,v);
psi=psi+Omega*At;
x=[0,psi];
o
RES
NMO
Total_Power=sum(RES.Power);
Total_Thrust=sum(RES.Thrust);
Total_Torque=sum(RES.Torque);
Power_vectorY(o)=Total_Power;
Power_thrustY(o)=Total_Thrust;
Power_TorqueY(o)=Total_Torque;
t=t+At;
WT.tt=0.4*sin(10*t);
end
plot(Time,Power_thrustY)