clear all
WT.nB=3;
WT.R=0.61;
xs=0.08;
xt=WT.R;
rstep=0.04;
WT.r=[xs:rstep:xt];
WT.chord=[0.1423 0.127 0.106 0.0875 0.0715 0.0575 0.0466 0.0402 0.037 0.035 0.0327 0.0304 0.028 0.02571]; %chord (m)
%WT.twist=[9.85 11.5 14.6 16.4 18.8 20.2 20.2 20.2 21.7 23.6 25 25.8 26.5 27.03]*pi/180; %twist angle (rad)
%WT.twist=[17.18 15.53 12.43 10.63 8.23 6.83 6.83 6.8 5.33 3.43 2.03 1.23 0.53 0]*pi/180;
WT.twist=[27.03 28.68 27.03 23.93 22.13 19.73 18.33 18.33 18 16 14.93 13.53 12.73 12.03 11.5]*pi/180;
%Wind conditions 
Sim.rho=1.225 ;
Sim.KinVisc=1.470e-5 ; % Kinematic viscosity [m^2/s] (for Reynolds number)
Sim.PITCH=0.0 ; % Pitch angle [rad] |+

Sim.RPM=650 ; % Rotational velocity [rad/s]
Wind.V0=12.5;% wind velocity
%algorithm parameters 
Algo.nbIt=60 ; % Maximum number of iterations
Algo.aTol= 0.01; % Tolerance in axial induction
Algo.bTipLoss=true; % True if tip -losses are applied

fBEMsteady(WT,Sim,Wind,Algo)
