
function [RES] = fBEMsteady(WT ,Sim ,Wind , Algo)
% fBEMsteady: steady BEM code implementation
% Author: E. Branlard
% --- Necessary parameters
% Environment and operating conditions
num=0;
V0 = Wind.V0 ; % Incoming Wind [m/s]
rho = Sim.rho ; % Air density [kg/m^3]
KinVisc = Sim.KinVisc ; % Kinematic viscosity [m^2/s] (for Reynolds number)
pitch = Sim.PITCH ; % Pitch angle [rad]
Omega = Sim.RPM*2* pi /60 ; % Rotational velocity [rad/s]
% Rotor geometry
nB = WT.nB ; % Number of blades
R = WT.R ; % Rotor radius [m]
r = WT.r ; % Radial positions [m] (1 x ne)
chord = WT.chord ; % Chord [m] (1 x ne)
twist = WT.twist ; % Twist [rad] (1 x ne)
% Algorithm options
nbIt = Algo.nbIt ; % Maximum number of iterations
aTol = Algo.aTol ; % Tolerance in axial induction
bTipLoss = Algo.bTipLoss; % True if tip -losses are applied
% --- Derived parameters
lambda_r = Omega*r/V0 ; % Local tip -speed ratio
sigma = chord*nB./(2* pi*r); % Solidity
% -- Initialize result variable
RES.Pn= zeros (1, length (r)); % normal force per length
RES.Pt= zeros (1, length (r)); % tangential force per length
RES.A =zeros (1, length (r));
RES.rey =zeros (1, length (r));
RES.CLac =zeros (1, length (r));
RES.CDac =zeros (1, length (r));
RES.aac=zeros (1, length (r));
RES.cnac=zeros (1, length (r));
RES.ctac=zeros (1, length (r));
% --- Loop on blade elements
for e=1: length (r)
% --- Step 0: initial guess
a = 0.3*0 ;
aprime = 0.01*0;
% --- Iteration loop
for i=1:nbIt
% --- Step 1: Wind Components
Ut = Omega *r(e) *(1+aprime) ;
Un = V0*(1-a) ;
Vrel_norm = sqrt (Un.^2+Ut.^2) ;
Re = Vrel_norm* chord(e)/KinVisc; % Reynolds number
%Re=80000;
if Re>250000
    Re=250000;
    fprintf('enter reynolds\n')
end
if Vrel_norm>100
   fprintf('NAN velocity\n')
end
% --- Step 2: Flow Angle
phi = atan2 (Un ,Ut); % [rad]
fprintf('-------------------------\n')
if( imag (phi)~=0); fprintf ('Algorithm failed: r =%.2f\n',r(e));
break ; 
end
% --- Tip loss
F = fTipLoss(nB ,r,R,phi ,bTipLoss,e);
% --- Step 3: Angle of attack
alpha=phi -(twist (e)+ pitch); % [rad]
% --- Step 4: Airfoil coefficients (and dynamic stall)
[Cl , Cd] = fAeroCoeff(alpha ,Re,e);
% --- Step 5: airfoil , load coeff and circulation
% Normal and tangential coefficients
cn = Cl.* cos(phi)+Cd.* sin(phi); % cnNoDrag = Cl.*cos(phi);
ct = Cl.* sin(phi)-Cd.* cos(phi); % ctNoDrag = Cl.*sin(phi);
% Local thrust and torque from BET
Ct=Vrel_norm^2/V0 ^2*sigma(e) *cn;
Cq=Vrel_norm^2/V0 ^2*sigma(e) *ct;
% Circulation for one blade
Gamma_B=0.5* Vrel_norm* chord(e)*Cl;
% --- Step 6: Induction Coefficients
% Storing last values
a_last = a ;
aprime_last = aprime;
[a,aprime] = fInductionCoefficients(a_last ,Ct ,Cq ,F,phi,sigma(e), lambda_r(e),cn) ;
% --- Convergence Criteria
if (i>3 && (abs(a-a_last)+abs(aprime-aprime_last))<aTol); break ;
end

end % iterative loop for one element
num=num+1;
fprintf ('section %.2f\n',num)
if(i==nbIt); fprintf ('Maximum iterations reached at r =%.2f\n',r(e)); 
end
% --- Step 8: Aerodynamic Forces per length ( WITH DRAG)
RES.Pn(e) = 0.5*rho* Vrel_norm .^2*chord(e).*cn;
RES.Pt(e) = 0.5*rho* Vrel_norm .^2*chord(e).*ct;
RES.rey(e) = Re;
RES.A(e)=alpha;
RES.CLac(e) =Cl;
RES.CDac(e) =Cd;
RES.cnac(e)=cn;
RES.ctac(e)=ct;
RES.aac(e)=a;
end % loop on blade elements
RES.Thrust = nB* trapz (r, RES.Pn) ;
RES.Power = nB* trapz (r,r.*RES.Pt)*Omega; % NOTE: Trapz not optimal!
RES.CP = RES.Power /(0.5*rho*V0^3* pi*R^2);
RES.CT = RES.Thrust /(0.5*rho*V0^2* pi*R^2);
end
function [F]= fTipLoss(nB ,r,R,phi , bTipLoss, e , varargin)
% - Compute tip -loss factor
% NOTE: Many implementations possible! Minimalistic example:
F=1;
if (bTipLoss & sin(phi) >0.01)
F=2/ pi* acos ( exp(-nB/2*(R-r(e))/(r(e)* sin(phi))));
end
end
function [Cl , Cd] = fAeroCoeff(alpha ,Re , e, varargin)
% - Interpolation of polars for a given alpha and Reynolds number
% - Dynamic stall implementationn

 alp=alpha*180/pi;
 alpha_st=10*pi/180;
 
 if alp>=10
    PN=false;
    fprintf ('Stall \n');
    [CLp,CDp,CLL,CDD]=airfoilinterpol(Re,PN,e);
    L1=length(CLL');
    L2=length(CDD');
    if(length(CLL)<8 || length(CDD)<8 || nnz(~CLL)>0 || nnz(~CDD)>0)
      [CLp,CDp,CLL,CDD]=airfoilinterpol(round(Re,-3),PN,e);
      fprintf ('Recalculation\n'); 
       if(length(CLL)<8 || length(CDD)<8 || nnz(~CLL)>0 || nnz(~CDD)>0)
        [CLp,CDp,CLL,CDD]=airfoilinterpol(round(Re,-4),PN,e);
        fprintf ('2dRecalculation\n'); 
       end

    end
    Clstall=polyval(CLp,10);
    Cdstall=polyval(CDp,10);
    Cl=sin(2*round(alpha,2))+(Clstall-sin(2*round(alpha_st,1)))*((sin(alpha_st))/(cos(round(alpha_st,1))^2))*((cos(round(alpha,2))^2)/(sin(round(alpha,2))));
    Cd=2*sin(round(alpha,2))^2+(Cdstall-2*sin(round(alpha_st,1))^2)*(cos(round(alpha,2))/cos(round(alpha_st,1)));

  elseif alp<=0
    fprintf ('Negative\n');
    PN=true;
    [CLp,CDp,CLL,CDD]=airfoilinterpol(Re,PN,e);
    if(length(CLL)<8 || length(CDD)<8 || nnz(~CLL)>0 || nnz(~CDD)>0)
      [CLp,CDp,CLL,CDD]=airfoilinterpol(round(Re,-3),PN,e);
      fprintf ('Recalculation\n'); 
      if(length(CLL)<8 || length(CDD)<8 || nnz(~CLL)>0 || nnz(~CDD)>0)
      [CLp,CDp,CLL,CDD]=airfoilinterpol(round(Re,-4),PN,e);
      fprintf ('2dRecalculation\n'); 
      end
    end
   Cl=polyval(CLp,round(alp,1));
   Cd=polyval(CDp,round(alp,1));
 else
    PN=false;
    [CLp,CDp,CLL,CDD]=airfoilinterpol(Re,PN,e);
    if(length(CLL)<8 || length(CDD)<8 || nnz(~CLL)>0 || nnz(~CDD)>0)
      [CLp,CDp,CLL,CDD]=airfoilinterpol(round(Re,-3),PN,e);
      fprintf ('Recalculation\n'); 
       if(length(CLL)<8 || length(CDD)<8 || nnz(~CLL)>0 || nnz(~CDD)>0)
      [CLp,CDp,CLL,CDD]=airfoilinterpol(round(Re,-4),PN,e);
      fprintf ('2dRecalculation\n'); 
      end
    end
 Cl=polyval(CLp,round(alp,1));
 Cd=polyval(CDp,round(alp,1));
 end
end

% function [Cl , Cd] = fAeroCoeff(alpha ,Re , varargin)
% Cl=0;
% Cd=0;
% % - Interpolation of polars for a given alpha and Reynolds number
% % - Dynamic stall implementation
%  alp=alpha*180/pi;
%  alpha_st=13*pi/180;
%  if alp>=13
%      [Clstall,Cdstall]=airfoilCpCd(Re,13);
%       if (isempty(Clstall))
%          fprintf ('Entered Clstall Cdstall recalculation\n')
%          [Clstall,Cdstall]=airfoilCpCd(round(Re,-3),13);
%          if (isempty(Clstall))
%             fprintf ('Entered 2nd Clstall Cdstall recalculation\n')
%             [Clstall,Cdstall]=airfoilCpCd(round(Re,-4),13);
%             if (isempty(Clstall))
%               fprintf ('Entered 3nd Clstall Cdstall recalculation\n')
%               [Clstall,Cdstall]=airfoilCpCd(round(Re,-5),13);
%             end
%          end
%      end
%      Cl=sin(2*round(alpha,2))+(Clstall-sin(2*round(alpha_st,1)))*((sin(alpha_st))/(cos(round(alpha_st,1))^2))*((cos(round(alpha,2))^2)/(sin(round(alpha,2))));
%      Cd=2*sin(round(alpha,2))^2+(Cdstall-2*sin(round(alpha_st,1))^2)*(cos(round(alpha,2))/cos(round(alpha_st,1)));
%  elseif alp<12
%      [Cl,Cd]=airfoilCpCd(Re,round(alp,2));
%      if (isempty(Cl))
%          fprintf ('Entered Cl Cd recalculation\n')
%          [Cl,Cd]=airfoilCpCd(round(Re,-3),round(alp,0));
%      end
%  end
% end
function [a,aprime] = fInductionCoefficients( a_last ,Ct ,Cq ,F,lambda_r,phi,sigma,cn , varargin) 
% - Compute a, aprime and the local thrust coefficient Ct
% - Perform High - thrust correction (e.g. a-Ct relation)
% - Perform relaxation on axial induction ( only if steady simulation)
% - Perform wake - rotation correction
% NOTE: Many implementations possible! Minimalistic example:
[a,Ct] = fCorrectionHighThrust (Ct ,F,phi,sigma,cn,varargin); % a-Ct relation
a = 0.3*a + (1-0.3)* a_last ; % Relaxation
aprime = Cq / (4*F*(1-a)*lambda_r) ; % tangential induction
end
function [a,Ct] = fCorrectionHighThrust (Ct ,F,phi,sigma,cn,varargin)
% - Returns a and Ct applying the High - Thrust correction corregir
% NOTE: Many implementations possible! Minimalistic example:
k = [0.00,0.251163,0.0544955,0.0892074];
Ctb = Ct./F;
a = k(4)*Ctb.^3+k(3)*Ctb.^2+k(2)*Ctb+k(1);
ac=0.29;
K=(4*F*(sin(phi))^2)/(sigma*cn);
% apru=1/2*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2)^2+4*(K*ac^2-1)));
end