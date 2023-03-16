
function [RES,NMO] = fBEMUnsteady(WT ,Sim ,Wind , Algo,NMO,x,v,p)
%NMO (n-1) in time step
% Rotor geometry
nB = WT.nB ; % Number of blades
R = WT.R ; % Rotor radius [m]
r = WT.r ; % Radial positions [m] (1 x ne)
rhub=WT.rhub;
chord = WT.chord ; % Chord [m] (1 x ne)
twist = WT.twist ; % Twist [rad] (1 x ne)
dr=WT.dr;
fs=0;


%wind and  Sim characteristic 
V0 = Wind.V0 ; % Incoming Wind [m/s]
rho = Sim.rho ; % Air density [kg/m^3]
KinVisc = Sim.KinVisc ; % Kinematic viscosity [m^2/s] (for Reynolds number)
pitch = Sim.PITCH ; % Pitch angle [rad]
Omega = Sim.RPM*2* pi /60 ; % Rotational velocity [rad/s]
dt=Sim.dt;

bTipLoss = Algo.bTipLoss; % True if tip -losses are applied


%%--initialization
n=[0;0;-1]; %thrust direction 
n_rotor=[0;0;1];%normal direction to the rotor 
ne=length(WT.r);%number of elements 
W0=ones(3,ne,nB);%induced velocity without yaw an tilt
W=zeros(3,ne,nB);%induced velocity final 
W_qs=zeros(3,ne,nB);%induced velocity quasistatic
W_int=zeros(3,ne,nB);%induced velocity intermediate 


% -- Initialize result variable
RES.Pn= zeros (nB, length (r)); % normal force per length
RES.Pt= zeros (nB, length (r)); % tangential force per length
RES.A =zeros (nB, length (r));%angle of attack per airfoil
RES.rey =zeros (nB, length (r));%reynold per airfoil
RES.CLac =zeros (nB, length (r));%Cl per airfoil
RES.CDac =zeros (nB, length (r));%Cd per 
RES.L =zeros (nB, length (r));%Cl per airfoil
RES.D =zeros (nB, length (r));%Cd per airfoil
RES.aac=zeros (nB, length (r));%a factor acumulation
RES.cnac=zeros (nB, length (r));%normal coefficient 
RES.ctac=zeros (nB, length (r));%tangential coefficient
RES.BladeTorque=zeros(1,nB);%blade torque
RES.BladeThrust=zeros(1,nB);%blade thrust
RES.Thrust = zeros(1,nB);
RES.Power = zeros(1,nB);
RES.Torque = zeros(1,nB);
RES.CP = zeros(1,nB);
RES.CT = zeros(1,nB);
Pz=0;
Py=0;
khi=NMO.khi;%khi initalitation
psi0=0;%potition


ty=WT.ty;
tt=WT.tt;
tw=WT.tw;
tc=WT.tc;
Yo=[0;p;0];
RR=[cos(tt);cos((pi/2));-cos((pi/2)-tt);]*4;
Rn=cross(Yo,RR);

% Angles between blades updated depending on Dynamic degree of freedom 2
Vpsi0=mod(0:(2*pi/nB):(2*pi/nB)*(nB-1),2*pi);
Vpsi= mod(Vpsi0 + x(2),2*pi); % 


%rotation matrix
a12=[1 0 0;0 1 0;0 0 1]*[cos(tt) 0 -sin(tt); 0 1 0; sin(tt) 0 cos(tt)]*[1 0 0; 0 cos(ty) sin(ty); 0 -sin(ty) cos(ty)];%rotation matrix base to top
a23=[cos(tw) sin(tw) 0; -sin(tw) cos(tw) 0; 0 0 1];%rotation matrix top to hub 
a34=[cos(tc) 0 -sin(tc); 0 1 0; sin(tc) 0 cos(tc)];%rotation matrix hub to blade_point
a14=a12*a23*a34;%rotation matrix base to blade_point
a13=a12*a23; %rotation matrix base to hub
a41=a12'*a23'*a34';%rotation matrix blade_point to base
a31=a12'*a23';%rotation matrix top to blade_point
a21=a12';%rotation matrix hub to blade_point



% --- Derived parameters

sigma = chord*nB./(2* pi*r); % Solidity



%----loop on number of blades
for l=1:nB
psi=Vpsi(l)
a23=[cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
a14=a12*a23*a34;%rotation matrix base to blade_point
a41=a14';%rotation matrix blade_point to base
a31=a12'*a23';%rotation matrix top to blade_point
a13=a31'; %rotation matrix base to hub
a21=a12';%rotation matrix hub to blade_point

% --- Loop on blade elements
for e=1: ne


    % --- Step 0: relative wind
    rb_in4=[r(e); 0; 0];
    rb_in1=a41*rb_in4;
    rs_in1=a21*WT.rs_in1;
    r_position=WT.rt_in1+rs_in1+rb_in1;
    V0_in1=[0,0,V0]'-Rn;%incoming wind seen by the hole system
    V0_in4=a14*V0_in1;%incoming wind seen by the blade
    V0_in3=a13*V0_in1;%incoming wind seen by the hub 
    V0_in4=[0;V0_in4(2);V0_in4(3)];% x component neglected
    V0_in3=[0;V0_in3(2);V0_in3(3)];% x component neglected
    Vb_in4=[0;  -Omega*r(e); 0]; %blade radial speed for r 
    Vrel_in4=V0_in4+NMO.W(:,e,l)+Vb_in4;%relative wind velocity 
    lambda_r = Omega*r(e)/norm(V0_in3) ; % Local tip -speed ratio
    %----Step 1: Wind components
    Ut = -Vrel_in4(2);
    Un = Vrel_in4(3) ;
    Vrel_norm = sqrt (Un.^2+Ut.^2) ;
    Re = Vrel_norm* chord(e)/KinVisc; % Reynolds number
    if Re>250000
        Re=250000;
        fprintf('enter reynolds\n')
    end
    if Vrel_norm>100
        fprintf('NAN velocity\n')
    end


    % --- Step 2: Flow Angle
    phi = atan2 (Un ,Ut); % [rad] flow angle
    
    % --- Tip loss
    F = fTipLoss(nB ,r,R,phi ,bTipLoss,e,rhub);
    
    
    % --- Step 3: Angle of attack
    alpha=phi -(twist (e)+ pitch); % [rad]
    
    
    % --- Step 4: Airfoil coefficients (and dynamic stall)
    [Cl , Cd] = fAeroCoeff(alpha ,Re,e);


    % --- Step 5: airfoil , load coeff and circulation
    % Normal and tangential coefficients
    cn = Cl.* cos(phi)+Cd.* sin(phi); % cnNoDrag = Cl.*cos(phi);
    ct = Cl.* sin(phi)-Cd.* cos(phi); % ctNoDrag = Cl.*sin(phi);
    Ct=Vrel_norm^2/V0 ^2*sigma(e) *cn;
    Cq=Vrel_norm^2/V0 ^2*sigma(e) *ct;
    % Local thrust and torque from BET
     %     Ct=Vrel_norm^2/V0 ^2*sigma(e) *cn;
     %     Cq=Vrel_norm^2/V0 ^2*sigma(e) *ct;
    % Circulation for one blade
    Gamma_B=0.5* Vrel_norm* chord(e)*Cl;

    L=0.5*rho*norm(Vrel_in4).^2*chord(e)*Cl;

    %-----induction factors on coordinate system
    Wn_in4=[0 ; 0 ; NMO.W(3,e,l) ];
    Wn_in3= a34'*Wn_in4;
    nnW_in3= n.*(n .*Wn_in3); 
    V_prime_induction_in3= V0_in3+nnW_in3;
    sign=1;
    if(V_prime_induction_in3(3)<0)
        sign=-1;
    end
    a=(norm(V0_in3)-sign*norm(V_prime_induction_in3))/norm(V0_in3);%induction factor for zero yaw 
    a_last=a;
    [a,aprime] = fInductionCoefficients(a_last ,cn ,ct ,F,phi,sigma(e), lambda_r,Vrel_norm,V0_in4,Wn_in4,V0_in3) ;
    W_y_qs = - Omega*r(e)*aprime;
    W_z_qs = -norm(V0_in3)*a;
    W_qs(:,e,l)=[0 ; W_y_qs; W_z_qs];


    %dynamic wake 
    tau1=(1.1*R)/((1-1.3*min(a,0.5))*norm(V0_in4));%time constants Oye model
    tau2=(0.39-0.26*((r(e)/R)^2))*tau1;%time constants
    H=W_qs(:,e,l)+0.6*tau1*((W_qs(:,e,l)-NMO.W_qs(:,e,l))/dt);
    W_int(:,e,l)=H+(NMO.W_int(:,e,l)-H)*exp(-dt/tau1);
    W0(:,e,l)=W_int(:,e,l)+(NMO.W(:,e,l)-W_int(:,e,l))*exp(-dt/tau2);
    %Yaw/tilt model
%     W(:,e,l)=W0(:,e,l);%non yaw and tilt 
     V0_in2=a12*V0_in1;
     psi0=atan2(V0_in2(2),V0_in2(1));
     meanW_in4=[0;0;mean(W0(3,e,:))];
     meanW_in2=a34'*meanW_in4;
     V_primeK=V0_in2+meanW_in2;
     khi=acos(dot(n_rotor,V_primeK)/norm(V_primeK));
     W(:,e,l)=W0(:,e,l)*(1+(r(e)/R)*tan(khi/2)*cos(Vpsi(l)-psi0));

% --- Step 8: Aerodynamic Forces per length ( WITH DRAG)
    RES.Pn(l,e) = 0.5*rho* Vrel_norm .^2*chord(e).*cn;
    RES.Pt(l,e) = 0.5*rho* Vrel_norm .^2*chord(e).*ct;
    RES.L(l,e)=0.5*rho*Vrel_norm.^2*chord(e)*Cl;
    RES.D(l,e)=0.5*rho*Vrel_norm.^2*chord(e)*Cd;
    RES.rey(l,e) = Re;
    RES.A(l,e)=alpha;
    RES.CLac(l,e) =Cl;
    RES.CDac(l,e) =Cd;
    RES.cnac(l,e)=cn;
    RES.ctac(l,e)=ct;
    RES.aac(l,e)=a;
   
end % loop on blade elements   


fprintf('-------------------------\n')
RES.Thrust(l) = trapz (r, RES.Pn(l,:)) ;
RES.Power(l) = trapz (r,r.*RES.Pt(l,:))*Omega;
RES.Torque(l) = trapz (r,r.*RES.Pt(l,:));% NOTE: Trapz not optimal!
% RES.CP(l) = RES.Power /(0.5*rho*V0^3* pi*R^2);
% RES.CT (l)= RES.Thrust /(0.5*rho*V0^2* pi*R^2);




end%loop on blades
NMO.W_qs=W_qs;
NMO.W0=W0;
NMO.W=W;
NMO.W_int=W_int;
NMO.khi=khi;
NMO.fs=fs;
end
function [F]= fTipLoss(nB ,r,R,phi , bTipLoss, e ,rhub, varargin)
% - Compute tip -loss factor
% NOTE: Many implementations possible! Minimalistic example:
F_t=1;
if (bTipLoss && sin(phi) >0.01)
F_t=2/ pi* acos ( exp(-nB/2*(R-r(e))/(r(e)* sin(phi))));
end
 F_n=2/pi*acos(exp(-nB/2*(r(e)-rhub)/(rhub*sin(phi))));
 F=F_t*F_n;
end
function [Cl , Cd] = fAeroCoeff(alpha ,Re , e, varargin)
% - Interpolation of polars for a given alpha and Reynolds number
% - Dynamic stall implementationn

 alp=alpha*180/pi;
 alpha_st=10*pi/180;
 
 if alp>=10
    PN=false;
    %%fprintf ('Stall \n');
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


function [a,aprime] = fInductionCoefficients(a_last ,cn ,ct ,F,phi,sigma, lambda_r,Vrel,V0_in4,Wn_in4,V0_in3, varargin) 
ac=0.2;                                       
if a_last<=ac
    fg=1;
elseif a_last>ac
    fg=ac/a_last*(2-ac/a_last);
end
%     W_z_qs=-nB*BEM.L*cosd(phi)/(4*pi*rho*r(e)*BEM.F*norm(V0_in4+fg*nnW_in4 ));
%     W_y_qs=-nB*BEM.L*sind(phi)/(4*pi*rho*r(e)*BEM.F*norm(V0_in4+fg*nnW_in4 ));
    W_z_qs=Vrel.^2.*cn/(4*F*norm(V0_in4+fg*Wn_in4 ) )*sigma;
    W_y_qs=Vrel.^2.*ct/(4*F*norm(V0_in4+fg*Wn_in4 ) )*sigma;
    a=W_z_qs/norm(V0_in3);
    aprime=W_y_qs/(lambda_r*norm(V0_in3));
end
