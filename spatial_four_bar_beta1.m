% SPATIAL  FOUR-BAR  MECHANISM OF RSSR TYPE WITH DYNAMIC FRICTION. POSITION 
% AND FORCE ANALYSIS
%
% Version: Beta1
%
% Matlab  algorithm implementing the analytic solutions of the RSSR linkage 
% with dynamic Coulomb friction. Based on paper:
%
% Cremasco  Coelho,  G.,  “Exact  solutions for the spatial  four-bar  with 
% dynamic  friction:  Position,  force,  Painlevé paradox,  and singularity 
% asymptotics”, Mech. Mach. Theory 142 (2025), 106163.
%
% Implemented in Matlab R2018b.
%
% Remark: 
%
% - The  primary  goal  of this script was to produce numeric evidence that
%   the solution on the paper is correct.  It  was  used to produce many of
%   the figures shown there.
% - The current algorithm was tested for some RSSR configurations. Particu-
%   larly,  the inputs originally adopted here will reproduce results shown
%   in  the  paper,  which  have  shown  agreement  with  ADAMS   multibody
%   simulations. Nonetheless, NOT  ALL  RSSR  CONFIGURATIONS WERE TESTED in 
%   this computer implementation.  In  other  words, more testing should be
%   done to guarantee it works as a general-purpose RSSR solver.
% - If  you  identify  errors,  please  notify the author. The email can be 
%   found on the paper.
% - A huge part of this m-file contains raw data from ADAMS simulations for
%   the  RSSR  mechanism simulated in the article. This was used to compare
%   exact  and  ADAMS  solutions.  Unless  you want these results, for some 
%   reason, you can delete all the final lines indicated by the end of this
%   script.
%
% Guilherme Cremasco Coelho
% sept-10th-2026

%--------------------------------------------------------------------------
% START
% INPUTS
%--------------------------------------------------------------------------
% Arm lengths [mm]
A = 201.6318;
B = 350.0;
C = 201.6318;
D = 403.1627;

% Inclination of body 3 (C) [deg]
theta = 7.294;
theta=theta*pi/180;

% Initial position of vector A [mm]
vec_A0=A.*[-0.379646;0.400777;0.833815];

% Vector D [mm]
vec_D=D.*[0.524328;0.835662;-0.163552];

% Unit vector e1' [-]
vec_e1_prim=[0.009243;-0.951003;0.309043];

% Input Angle [deg]
alpha_a=-0.0:-1:-236;
alpha_a=alpha_a*pi/180;

% Select the solution branch (1 or 2)
branch=2;

% There exists only two types of boundary conditions.
% BC I: load applied on joint I and velocity on joint J
% BC J: load applied on joint J and velocity on joint I
% Only signs are specified, for load or for velocity:
% S_M_INPUT: sign for torque M from ground to bellcrank in joint I
% S_T_INPUT: sign for torque T from ground to bellcrank in joint J
% S_I: sign of angle alpha_a (joint I)
% S_J: sign of angle alpha_c (joint J)
BC_TYPE='J';S_T_INPUT=+1;S_I=+1;
%BC_TYPE='I';S_M_INPUT=+1;S_J=+1; 

% Friction-related data for revolute i
% Remark:  the Revolute Joint Friction Model adopted here is NOT  THE  SAME
% as those implemented in ADAMS. Details can be found in the paper.
R_i=5;       % [mm]
R_I=5;       % [mm]
R_I_ast=5;   % [mm]
mu_i=.7;     % [-]
mu_I=.7;     % [-]
mu_I_ast=.7; % [-]
h_i=30;      % [mm]

% Friction-related data for revolute j
R_j=5;       % [mm]
R_J=5;       % [mm]
R_J_ast=5;   % [mm]
mu_j=.6;     % [-]
mu_J=.6;     % [-]
mu_J_ast=.6; % [-]
h_j=30;      % [mm]

% Friction-related data for spherical joint m
R_m=20;      % [mm]
mu_m=0.2;    % [-]

% Friction-related data for spherical joint n
R_n=20;      % [mm]
mu_n=0.8;    % [-]
%--------------------------------------------------------------------------
% END
% INPUTS
%--------------------------------------------------------------------------


% Usual euclidean basis vectors e1, e2, e3
vec_e1=[1;0;0];
vec_e2=[0;1;0];
vec_e3=[0;0;1];

% Fix eventual normalization issues for e1'
vec_e1_prim=vec_e1_prim./norm(vec_e1_prim);
n1=vec_e1_prim(1);n2=vec_e1_prim(2);n3=vec_e1_prim(3);

% Compute coefficients.
% Eqs (30), (37), (43) 
c_i=R_i*mu_i;
c_I=R_I*mu_I;
c_I_ast=R_I_ast*mu_I_ast;
c_j=R_j*mu_j;
c_J=R_J*mu_J;
c_J_ast=R_J_ast*mu_J_ast;
c_m=R_m*mu_m;
c_n=R_n*mu_n;

% Starting  the  variable  that  stores the number of solutions for a given
% mechanism position
N_SOLUTIONS=0;

% Running over angle alpha_a (the angle of body 1 <A>):
for k=1:1:length(alpha_a)
    ca=cos(alpha_a(k));sa=sin(alpha_a(k));
    
    % Computes rotation matrix R
    % Eq (1)
    R11=ca+n1^2*(1-ca);
    R12=n1*n2*(1-ca)-n3*sa;
    R13=n2*sa+n1*n3*(1-ca);
    R21=n3*sa+n1*n2*(1-ca);
    R22=ca+n2^2*(1-ca);
    R23=-n1*sa+n2*n3*(1-ca);
    R31=-n2*sa+n1*n3*(1-ca);
    R32=n1*sa+n2*n3*(1-ca);
    R33=ca+n3^2*(1-ca);
    R=[R11,R12,R13;R21,R22,R23;R31,R32,R33];
    
    % Computes the vector A=R*A0
    % Eq (1)
    aux=R*vec_A0;
    vec_A(1:3,k)=aux;
    
    % Coefficients of the characteristic polynomial:
    % Eq (6)
    k2=(aux(1)-vec_D(1))^2+(aux(2)-vec_D(2))^2+(aux(3)-vec_D(3))^2-B^2+C^2+...
        2*C*((aux(1)-vec_D(1))*sin(theta)+(aux(2)-vec_D(2))*cos(theta));
    k1=-4*C*cos(theta)*(aux(3)-vec_D(3));
    k0=(aux(1)-vec_D(1))^2+(aux(2)-vec_D(2))^2+(aux(3)-vec_D(3))^2-B^2+C^2+...
        2*C*((aux(1)-vec_D(1))*sin(theta)-(aux(2)-vec_D(2))*cos(theta));
    
    % Quadratic determinant
    delta_2=k1^2-4*k0*k2;
    
    % Performs basic assembly analysis based on determinant and 
    % coefficients. Also computes the two mechanism solutions.
    % Eq (5)
    if k2~=0 & delta_2>0
        tau_1=(-k1+sqrt(delta_2))/(2*k2);
        tau_2=(-k1-sqrt(delta_2))/(2*k2);
    elseif k2~=0 & delta_2==0
        tau_1=-k1/(2*k2);
        tau_2=tau_1;
    elseif k2~=0 & k1~=0
        tau_1=-k0/k1;
        tau_2=tau_1;
    elseif k2~=0 & delta_2<0
        disp(['The mechanism is not solvable for k=',num2str(k)]);return;
    elseif k2==0 & k1==0 & k0~=0
        disp(['The mechanism is not solvable for k=',num2str(k)]);return;
    elseif k2==0 && k1==0 & k0==0
        disp(['The mechanism is degenerate for k=',num2str(k)]);return;
    end
    
    % Selects the mechanism branch according to the specified by the user.
    if branch==1
        tau=tau_1;
    else
        tau=tau_2;
    end
    
    % Recovers the body 3 (C) output angle from its tan-half-angle tau.
    alpha_c(k)=2*atan(tau);

    % Reconstructs vector C
    % Eq (3)
    aux=C.*[-sin(theta);...
        cos(theta).*(1-tau^2)/(1+tau^2);...
        cos(theta).*(2*tau)/(1+tau^2)];
    vec_C(1:3,k)=aux;

    % Reconstructs vector B from loop equation B=C+D-A
    vec_B(1:3,k)=vec_C(1:3,k)+vec_D-vec_A(1:3,k);

    % Computes the frictionless gain T/M using the equation
    % T/M = [e1.(CxB)]/[e1'.(AxB)]
    % This equation is not in the paper, but can be obtained directly from
    % Eqs (14), (15) and (16)
    T_over_M(k)=-dot(vec_e1,cross(vec_C(1:3,k),vec_B(1:3,k)))/...
        dot(vec_e1_prim(1:3),cross(vec_A(1:3,k),vec_B(1:3,k)));
    
    % Straightforward way to get the relation d(alpha_a)/d(alpha_c):
    % Instead  of  using  numeric  derivatives  or  super-lengthy analytic
    % "Mathematica" expressions, the Virtual Work Principle can be invoked
    % to produce the compact relation below:
    % (The  equation  was  removed from the paper as a means to address  a
    % reviewer comment).
    da_dc=-T_over_M(k);
 
    if BC_TYPE=='I'
        % IF the given boundary condition are:
        % Torque in I (M) and Velocity in J (S_J):
        S_M=S_M_INPUT;
        S_I=sign(da_dc).*S_J;
        
        % S_T for the frictionless case:
        S_T=S_M.*sign(T_over_M(k));
    elseif BC_TYPE=='J'
        % IF the given boundary condition are:
        % Torque in J (T) and Velocity in I (S_I):
        S_T=S_T_INPUT;
        S_J=sign(da_dc).*S_I;

        % S_M for the frictionless case:
        S_M=S_T_INPUT.*sign(T_over_M(k));
    end
    
    % The  problem  is  quasi-static with dynamic friction.  Consequently,
    % the amplitude of angular  velocity  is not important, just its sign.
    % We choose abs(d(alpha_a)/dt)=1, wlog.
    d_alpha_a_dt=S_I;
    
    % Straightforward chain rule application
    % Eq (44)
    d_alpha_c_dt=d_alpha_a_dt/da_dc;
    
    % Computes vectors e2' and e3'
    % Eqs (31)
    vec_e2_prim=vec_A(1:3,k)-...
        dot(vec_e1_prim(1:3),vec_A(1:3,k)).*vec_e1_prim(1:3);
    vec_e2_prim=vec_e2_prim./norm(vec_e2_prim);
    vec_e3_prim=cross(vec_e1_prim,vec_e2_prim);

    % Computes vector PHI using Eq (D4).
    % IMPORTANT REMARK:
    % In the published pdf version, there is a typo: "C" was placed outside 
    % the fraction. The published online html version is correct.  The SSRN
    % preprint is also corect,  though  the  article structure and equation
    % numbering is very different).
    vec_PHI=C/B*d_alpha_c_dt.*vec_e1-A/B*d_alpha_a_dt.*vec_e1_prim;
    
    % Computes vectors e1'', e2'' and e3''
    % Eq (20)
    vec_e1_prim_prim=vec_B(1:3,k);
    vec_e1_prim_prim=vec_e1_prim_prim/norm(vec_e1_prim_prim);
    vec_e2_prim_prim=vec_PHI./norm(vec_PHI);
    vec_e3_prim_prim=cross(vec_e1_prim_prim,vec_e2_prim_prim);

    % Computes the frictionless ratio F/M
    F_over_M(k)=-1./dot(vec_e1_prim,cross(vec_A(1:3,k),vec_e1_prim_prim));
    
    % Computes parameters gamma_a, gamma_c, gamma'
    % Eq (D3)
    gamma_a=d_alpha_a_dt.*dot(vec_e1_prim,vec_e1_prim_prim);
    gamma_c=d_alpha_c_dt.*dot(vec_e1,vec_e1_prim_prim);
    gamma_prim=(gamma_a-gamma_c)/2;
    
    % Computes parameters r' and s'
    % Eq (D3)
    r_prim=(d_alpha_a_dt.*dot(vec_e1_prim,vec_e2_prim_prim))^2+...
        (norm(vec_PHI)-d_alpha_a_dt.*dot(vec_e1_prim,vec_e3_prim_prim))^2;
    r_prim=sqrt(r_prim);
    
    s_prim=(d_alpha_c_dt.*dot(vec_e1,vec_e2_prim_prim))^2+...
        (norm(vec_PHI)-d_alpha_c_dt.*dot(vec_e1,vec_e3_prim_prim))^2;
    s_prim=sqrt(s_prim);
    
    % When gamma' is not zero, we compute r and s
    if gamma_prim~=0
        r=r_prim/gamma_prim;
        s=s_prim/gamma_prim;
    end
    
    calculate_N=1;
    
    % BEGIN: spin (d_csi_dt) calculation.
    % -----------------------------------------------------
    % Identify which is the case:
    % - Cases 1 and 2 are more general.
    % - Cases cases A thru G are more particular. Perhaps it is possible to
    % show  that  they never occur for any RSSR,  but proving or disproving
    % that  was  beyond the scope of the paper.  For  this reason, they are
    % considered to be possible (though "rare").
    
    % CASE 1:
    if gamma_prim~=0 & c_m==c_n & c_m~=0 & r~=0 & s~=0
        % Eqs (A9)
        z=-(r+s)/(r-s);
        if abs(z) < 1
            z=z;
        else
            z=1./z;
        end
        
        % Eq (A5):
        d_csi_dt=gamma_prim*z+(gamma_a+gamma_c)/2;

    % CASE 2:
    elseif gamma_prim~=0 & c_m~=c_n & c_m~=0 & c_n~=0 & r~=0 & s~=0
        % Here,  determining  dcsi/dt  requires  solving  a quartic. In the 
        % paper  I've shown that the solution, however, is unique. Proof of 
        % the  correct  root selection criterion among the four ones can be 
        % seen on the paper. Here is the implementation.
        t1=(c_m^2*s^2-c_n^2*r^2)/(c_m^2-c_n^2);
        t2=(c_m^2*s^2+c_n^2*r^2)/(c_m^2-c_n^2);
        Z=roots([1,0,(t1-2),-2*t2,t1+1]);
        Q=c_m.*(Z-1)./sqrt((Z-1).^2+r.^2)+c_n.*(Z+1)./sqrt((Z+1).^2+s^2);
        z=0;
        for counter=1:length(Z)
            if abs(Q(counter))<1e-5 & abs(real(Z(counter)))<1
                z=real(Z(counter));
            end
            ALL_Z=real(Z).*(abs(Q)<0.1);    
        end
        d_csi_dt=gamma_prim*z+(gamma_a+gamma_c)/2;
        D_CSI_DT(1:4,k)=gamma_prim*ALL_Z+(gamma_a+gamma_c)/2;

    % Particular case (a)
    elseif c_m==0 & c_n==0
        disp('Case (a)');
        % d_csi_dt exists but is not unique 
        d_csi_dt=0;
        calculate_N=0;
        N=vec_e1_prim_prim;
        
    % Particular case (b)
    elseif gamma_prim==0 & (c_m+c_n)~=0
        disp('Case (b)');
        d_csi_dt=-(gamma_a+gamma_c)/2;
        calculate_N=0;
        N=vec_e1_prim_prim;
        
    % Particular case (c)
    elseif gamma_prim~=0 & c_m==0 & c_n~=0
        if s~=0
            disp('Case (c), s not zero.');
            d_csi_dt=gamma_prim-(gamma_a+gamma_c)/2;
        else
            disp('Case (c), s=0.');
            %Existence is not guaranteed
            d_csi_dt=gamma_prim-(gamma_a+gamma_c)/2;
        end
        
    % Particular case (d) 
    elseif gamma_prim~=0 & c_m~=0 & c_n==0
        if r~=0
            disp('Case (d), r not zero.');
            d_csi_dt=-gamma_prim-(gamma_a+gamma_c)/2;
        else
            disp('Case (d), r==0.');
            %Existence is not guaranteed
            d_csi_dt=-gamma_prim-(gamma_a+gamma_c)/2;
        end
        
    % Particular case (e) 
    elseif gamma_prim~=0 & c_m~=0 & c_n~=0 & r~=0 & s==0
        if c_m>c_n && r^2<=4((c_m/c_n)^2-1)
            disp('Case (e), solution exists.');
            d_csi_dt=gamma_prim*(-1+sqrt((r^2*c_n^2)/(c_m^2-c_n^2)))...
                -(gamma_a+gamma_c)/2;
        else
            disp('Case (e), solution does not exist.');
            %A solution does not exist
            break;
        end
    
    % Particular case (f) 
    elseif gamma_prim~=0 & c_m~=0 & c_n~=0 & r==0 & s~=0
        if c_m<c_n & s^2<=4((c_n/c_m)^2-1)
            disp('Case (f), solution exists.');
            d_csi_dt=gamma_prim*(1-sqrt((s^2*c_m^2)/(c_n^2-c_m^2)))...
                -(gamma_a+gamma_c)/2;
        else
            disp('Case (f), solution does not exist.');
            %A solution does not exist
            break;
        end
        
    % Particular case (g)
    elseif gamma_prim~=0 & c_m~=0 & c_n~=0 & r==0 & s==0
        if c_m==c_n
            disp('Case (f), solution exists, but is non-unique.');
            % Solution is non-unique
            d_csi_dt=-(gamma_a+gamma_c)/2;
        else
            disp('Case (f), solution does not exist.');
            %A solution does not exist
            break;
        end
    end
    % END: spin (d_csi_dt) calculation.
    % -----------------------------------------------------
    
    % Computes the loads for the case with friction.
    % -----------------------------------------------------
    if calculate_N==1
        % Computes the angular velocity vectors omega_12'' and omega_32''
        % Eqs (28)
        omega_12_1_prim_prim=...
            -d_csi_dt+d_alpha_a_dt*dot(vec_e1_prim,vec_e1_prim_prim);
        omega_12_2_prim_prim=...
            d_alpha_a_dt*dot(vec_e1_prim,vec_e2_prim_prim);
        omega_12_3_prim_prim=...
            -norm(vec_PHI)+d_alpha_a_dt*dot(vec_e1_prim,vec_e3_prim_prim);
        
        omega_32_1_prim_prim=...
            -d_csi_dt+d_alpha_c_dt*dot(vec_e1,vec_e1_prim_prim);
        omega_32_2_prim_prim=...
            +d_alpha_c_dt*dot(vec_e1,vec_e2_prim_prim);
        omega_32_3_prim_prim=...
            -norm(vec_PHI)+d_alpha_c_dt*dot(vec_e1,vec_e3_prim_prim);
        
        vec_omega_12=omega_12_1_prim_prim*vec_e1_prim_prim+...
            omega_12_2_prim_prim*vec_e2_prim_prim+...
            omega_12_3_prim_prim*vec_e3_prim_prim;
        
        vec_omega_21=-vec_omega_12;

        vec_omega_32=omega_32_1_prim_prim*vec_e1_prim_prim+...
            omega_32_2_prim_prim*vec_e2_prim_prim+...
            omega_32_3_prim_prim*vec_e3_prim_prim;

        vec_omega_23=-vec_omega_32;
        
        Recalculate_S_F=0;
        
        %-----------------------------------------------------------------
        % START
        % ATTEMPTS A FIRST SOLUTION TO THE SIGN PROBLEM
        %-----------------------------------------------------------------
        % The first guess is arbitrary. We pick S_F as the sign that F 
        % would have in the frictionless case.
        if BC_TYPE=='I'
            % IF the given boundary condition are:
            % Torque in I (M) and Velocity in J (S_J):
            S_M=S_M_INPUT;
            
            % Computes 
            % sign(d(alpha_a)/dt)=
            %   =  sign(d(alpha_a)/d(alpha_c)) * sign(d(alpha_c)/dt)
            S_I=sign(da_dc).*S_J;

            % Computes, for the frictionless case, sign(T):
            % sign(T)=sign(M)*sign(T/M)
            S_T=S_M.*sign(T_over_M(k));
        elseif BC_TYPE=='J'
            % IF the given boundary condition are:
            % Torque in J (T) and Velocity in I (S_I):
            S_T=S_T_INPUT;
            S_J=sign(da_dc).*S_I;

            % Computes, for the frictionless case, sign(M):
            S_M=S_T_INPUT.*sign(T_over_M(k));
        end

        % Determine the first guess of S_F 
        % Eq (14)
        F_frictionless(k)=-S_M*...
            dot(vec_e1_prim,cross(vec_A(1:3,k),vec_e1_prim_prim));
        S_F=sign(F_frictionless(k));

        % Computes the vector PSI
        % Eq (48)
        vec_PSI=-c_m*S_F*vec_omega_12./norm(vec_omega_12)-...
            c_n*S_F*vec_omega_32./norm(vec_omega_32);

        % Computes the vector N''
        % Eq (50)
        N_2_prim_prim=dot(vec_PSI,vec_e3_prim_prim)/B;
        N_3_prim_prim=-dot(vec_PSI,vec_e2_prim_prim)/B;
        N_1_prim_prim=sqrt(1-N_2_prim_prim^2-N_3_prim_prim^2);

        vec_N=N_1_prim_prim*vec_e1_prim_prim+...
            N_2_prim_prim*vec_e2_prim_prim+...
            N_3_prim_prim*vec_e3_prim_prim;

        % *****************************************
        % Start:
        % Computes many friction-related parameters:
        % *****************************************
        % Eqs (B4)
        D_Icom1=dot(vec_N,vec_e1_prim);
        C_icom2=dot(vec_N,vec_e2_prim);
        C_icom3=dot(vec_N,vec_e3_prim);
        D_Jcom1=dot(vec_N,vec_e1);
        C_jcom2=dot(vec_N,vec_e2);
        C_jcom3=dot(vec_N,vec_e3);
        
        % Eq (B6)
        aux=cross(vec_A(1:3,k),vec_N)...
            +c_m.*S_F.*vec_omega_21./norm(vec_omega_21);
        C_Icom2=1/h_i.*dot(aux,vec_e3_prim);
        C_Icom3=1/h_i.*dot(aux,vec_e2_prim);
        C_Icom1=dot(aux,vec_e1_prim);

        % Eqs (B6) (cont'd)
        aux=cross(vec_C(1:3,k),vec_N)...
            -c_n.*S_F.*vec_omega_23./norm(vec_omega_23);
        C_Jcom2=1/h_j.*dot(aux,vec_e3);
        C_Jcom3=1/h_j.*dot(aux,vec_e2);
        C_Jcom1=dot(aux,vec_e1);

        % Eqs (B8)
        K_I=sqrt(C_Icom2^2+C_Icom3^2);
        K_i=sqrt((C_icom2-C_Icom2)^2+(C_icom3+C_Icom3)^2);
        K_I_ast=abs(D_Icom1);

        % Eqs (B10)
        K_J=sqrt(C_Jcom2^2+C_Jcom3^2);
        K_j=sqrt((C_jcom2-C_Jcom2)^2+(C_jcom3+C_Jcom3)^2);
        K_J_ast=abs(D_Jcom1);
        % *****************************************
        % End:
        % Computes many friction-related parameters:
        % *****************************************

        S=sign(da_dc);
        S_i=S_I;
        S_j=S_J;

        % *****************************************
        % START:
        % Computes the resulting sign to see if it 
        % matches the first guess
        % *****************************************
        if BC_TYPE=='I'
            % IF the given boundary condition are:
            % Torque in I (M) and Velocity in J (S_J):
            % Eq (55)(a)
            den_M=S_i.*S_F*(c_I.*K_I+c_i.*K_i+c_I_ast*K_I_ast)-C_Icom1;
            F_over_M_friction_try1(k)=1./den_M;

            % Computed sign(F) 
            % (will be seen if it matches the first guess or not)
            new_S_F=sign(F_over_M_friction_try1(k)).*S_M_INPUT;

            S_M=S_M_INPUT;
            
            % Eq (55)(b)
            num_M=C_Jcom1+S_j*S_F*(c_J.*K_J+c_j.*K_j+c_J_ast*K_J_ast);
            T_over_M_friction_try1(k)=num_M./den_M;
            M_over_T_friction_try1(k)=...
                1./T_over_M_friction_try1(k);

            S_T=sign(T_over_M_friction_try1(k)).*S_M;

            % Eq (55)(c)
            F_over_T_friction_try1(k)=...
                F_over_M_friction_try1(k).*M_over_T_friction_try1(k);

        elseif BC_TYPE=='J'                
            % IF the given boundary condition are:
            % Torque in J (T) and Velocity in I (S_I):
            % Eq (55)(c)
            den_T=C_Jcom1+S_j*S_F*(c_J.*K_J+c_j.*K_j+c_J_ast*K_J_ast);
            F_over_T_friction_try1(k)=1./den_T;

            new_S_F=sign(F_over_T_friction_try1(k)).*S_T_INPUT;

            S_T=S_T_INPUT;

            % Eq (55)(b) inversed
            num_T=S_i.*S_F*(c_I.*K_I+c_i.*K_i+c_I_ast*K_I_ast)-C_Icom1;
            M_over_T_friction_try1(k)=num_T./den_T;
            S_M=sign(M_over_T_friction_try1(k)).*S_T;
            T_over_M_friction_try1(k)=...
                1./M_over_T_friction_try1(k);

            % Eq (55)(a)
            F_over_M_friction_try1(k)=...
                F_over_T_friction_try1(k).*T_over_M_friction_try1(k);
        end
        
        N_SOL_S_F=0;
        if new_S_F==S_F
            N_SOL_S_F=1;
        else
            % If the chosen sign this is not a solution
            F_over_T_friction_try1(k)=nan;
            F_over_M_friction_try1(k)=nan;
            M_over_T_friction_try1(k)=nan;
            T_over_M_friction_try1(k)=nan;
        end
        % *****************************************
        % END:
        % Computes the resulting sign to see if it 
        % matches the first guess
        % *****************************************
        %-----------------------------------------------------------------
        % END
        % ATTEMPTS A FIRST SOLUTION TO THE SIGN PROBLEM
        %-----------------------------------------------------------------
        
        
        %-----------------------------------------------------------------
        % START
        % ATTEMPTS A SECOND SOLUTION TO THE SIGN PROBLEM
        %-----------------------------------------------------------------
        % The second guess is opposite to the first guess. The calculations
        % are just a repetition of what was done before, but now using:
        S_F=-S_F;

        if BC_TYPE=='I'
            % IF the given boundary condition are:
            % Torque in I (M) and Velocity in J (S_J):
            S_M=S_M_INPUT;
            S_T=S_M.*sign(T_over_M);
        elseif BC_TYPE=='J'
            % If Actuator is placed in joint J, the torque T direction
            % is given and the torque M in joint I is sought
            S_T=S_T_INPUT;
            S_M=S_T_INPUT.*sign(T_over_M);
        end

        % Computes the vector PSI
        % Eq (48)
        vec_PSI=-c_m*S_F*vec_omega_12./norm(vec_omega_12)-...
            c_n*S_F*vec_omega_32./norm(vec_omega_32);

        % Computes the vector N''
        % Eq (50)
        N_2_prim_prim=dot(vec_PSI,vec_e3_prim_prim)/B;
        N_3_prim_prim=-dot(vec_PSI,vec_e2_prim_prim)/B;
        N_1_prim_prim=sqrt(1-N_2_prim_prim^2-N_3_prim_prim^2);

        vec_N=N_1_prim_prim*vec_e1_prim_prim+...
            N_2_prim_prim*vec_e2_prim_prim+...
            N_3_prim_prim*vec_e3_prim_prim;

        % *****************************************
        % Start:
        % Computes many friction-related parameters:
        % *****************************************
        % Eqs (B4)
        D_Icom1=dot(vec_N,vec_e1_prim);
        C_icom2=dot(vec_N,vec_e2_prim);
        C_icom3=dot(vec_N,vec_e3_prim);
        D_Jcom1=dot(vec_N,vec_e1);
        C_jcom2=dot(vec_N,vec_e2);
        C_jcom3=dot(vec_N,vec_e3);

        % Eq (B6)
        aux=cross(vec_A(1:3,k),vec_N)...
            +c_m.*S_F.*vec_omega_21./norm(vec_omega_21);
        C_Icom2=1/h_i.*dot(aux,vec_e3_prim);
        C_Icom3=1/h_i.*dot(aux,vec_e2_prim);
        C_Icom1=dot(aux,vec_e1_prim);

        % Eqs (B6) (cont'd)
        aux=cross(vec_C(1:3,k),vec_N)...
            -c_n.*S_F.*vec_omega_23./norm(vec_omega_23);
        C_Jcom2=1/h_j.*dot(aux,vec_e3);
        C_Jcom3=1/h_j.*dot(aux,vec_e2);
        C_Jcom1=dot(aux,vec_e1);

        % Eqs (B8)
        K_I=sqrt(C_Icom2^2+C_Icom3^2);
        K_i=sqrt((C_icom2-C_Icom2)^2+(C_icom3+C_Icom3)^2);
        K_I_ast=abs(D_Icom1);

        % Eqs (B10)
        K_J=sqrt(C_Jcom2^2+C_Jcom3^2);
        K_j=sqrt((C_jcom2-C_Jcom2)^2+(C_jcom3+C_Jcom3)^2);
        K_J_ast=abs(D_Jcom1);
        % *****************************************
        % End:
        % Computes many friction-related parameters
        % *****************************************

        
        % *****************************************
        % START:
        % Computes the resulting sign to see if it 
        % matches the first guess
        % *****************************************
        if BC_TYPE=='I'
            % IF the given boundary condition are:
            % Torque in I (M) and Velocity in J (S_J):
            % Eq (55)(a)
            den_M=S_i.*S_F*(c_I.*K_I+c_i.*K_i+c_I_ast*K_I_ast)-C_Icom1;
            F_over_M_friction_try2(k)=1./den_M;

            % Computed sign(F) 
            % (will be seen if it matches the second guess or not)
            new_S_F=sign(F_over_M_friction_try2(k)).*S_M_INPUT;

            S_M=S_M_INPUT;

            % Eq (55)(b)
            num_M=C_Jcom1+S_j*S_F*(c_J.*K_J+c_j.*K_j+c_J_ast*K_J_ast);
            T_over_M_friction_try2(k)=num_M./den_M;
            M_over_T_friction_try2(k)=...
                1./T_over_M_friction_try2(k);

            S_T=sign(T_over_M_friction_try2(k)).*S_M;

            % Eq (55)(c)
            F_over_T_friction_try2(k)=...
                F_over_M_friction_try2(k).*M_over_T_friction_try2(k);

        elseif BC_TYPE=='J'                
            % IF the given boundary condition are:
            % Torque in J (T) and Velocity in I (S_I):
            % Eq (55)(c)
            den_T=C_Jcom1+S_j*S_F*(c_J.*K_J+c_j.*K_j+c_J_ast*K_J_ast);
            F_over_T_friction_try2(k)=1./den_T;

            new_S_F=sign(F_over_T_friction_try2(k)).*S_T_INPUT;

            S_T=S_T_INPUT;

            % Eq (55)(b) inversed
            num_T=S_i.*S_F*(c_I.*K_I+c_i.*K_i+c_I_ast*K_I_ast)-C_Icom1;
            M_over_T_friction_try2(k)=num_T./den_T;
            T_over_M_friction_try2(k)=...
                1./M_over_T_friction_try2(k);

            S_M=sign(M_over_T_friction_try2(k)).*S_T;

            % Eq (55)(a)
            F_over_M_friction_try2(k)=...
                F_over_T_friction_try2(k).*T_over_M_friction_try2(k);
        end

        if new_S_F==S_F
            N_SOL_S_F=N_SOL_S_F+1;
        else
            F_over_T_friction_try2(k)=nan;
            F_over_M_friction_try2(k)=nan;
            M_over_T_friction_try2(k)=nan;
            T_over_M_friction_try2(k)=nan;
        end
        
        if N_SOL_S_F==0
            disp('-----------------------------------------------------------------------------------------------');
            disp('Painlevé Paradox.');
            disp('No Solution to the sign problem. Check if static equilibrium is possible or mechanism is jammed');
            disp(['alpha_a = ', num2str(alpha_a(k)*180/pi),' deg']);
            disp('-----------------------------------------------------------------------------------------------');
            N_SOLUTIONS(k)=0;
        elseif N_SOL_S_F==1
            N_SOLUTIONS(k)=1;
        elseif N_SOL_S_F==2
            disp('-----------------------------------------------------------------------------------------------');
            disp('Painlevé Paradox.');
            disp('Two solutions to the sign problem. Check if there exists a friction-induced instability. ');
            disp(['alpha_a = ', num2str(alpha_a(k)*180/pi),' deg']);
            disp('-----------------------------------------------------------------------------------------------');
            N_SOLUTIONS(k)=2;
        end
        % *****************************************
        % END:
        % Computes the resulting sign to see if it 
        % matches the second guess
        % *****************************************
        %-----------------------------------------------------------------
        % END
        % ATTEMPTS A SECOND SOLUTION TO THE SIGN PROBLEM
        %-----------------------------------------------------------------
    end
    
    % Stores the calculated d(csi)/dt (the body 2 spin)
    list_d_csi_dt(k)=d_csi_dt;

end


%--------------------------------------------------------------------------
% START
% PLOTTING THE EXACT RESULTS
%--------------------------------------------------------------------------
% 1. Plots the input-output angles
figure;
plot(alpha_a*180/pi,unwrap(alpha_c)*180/pi,'k-'); 
P1=gca;
xlabel('\alpha_c [deg]');
ylabel('\alpha_a [deg]');
title('Input-Output angles');

% 2. Plots normalized spin
figure;
plot(alpha_a(1:end-5)*180/pi,list_d_csi_dt(1:end-5)/S_I,'k-');
P2=gca;
xlabel('\alpha_a [deg]');
ylabel('d\xi/d\alpha_a [-]');
title('Normalized Spin d\xi/d\alpha_a');

% 3. Plots pseudo-roll angle
    % 3.a Calculates the PSEUDO-roll angle:
    csi_exact=cumsum(list_d_csi_dt/S_I.*gradient(alpha_a));

    % 3.b. Applies a certain initial condition
    % (constant of integration - arbitrary)
    csi_exact=csi_exact-csi_exact(1);
    CONST=0.12;
    csi_exact=csi_exact-CONST;

    % 3.c Plots pseudo-roll angle
    figure;
    plot(alpha_a(1:end-5)*180/pi,csi_exact(1:end-5)*180/pi,'k-');
    P3=gca;
    xlabel('\alpha_a [deg]');
    ylabel('\xi [deg]');
    title('Pseudo-roll angle \xi');

% 4. Plots M/T for the cases with and without friction
figure;
plot(alpha_a*180/pi,1./T_over_M,'k-');hold on;
plot(alpha_a(1:end-5)*180/pi,1./T_over_M_friction_try1(1:end-5),'r-');
hold on;
plot(alpha_a(1:end-5)*180/pi,1./T_over_M_friction_try2(1:end-5),'k-');
hold on;
P4=gca;
xlabel('\alpha_a [deg]');
ylabel('M/T [-]');
title({'M/T: Ratio between net input and net output torques.';...
    'Exact vs ADAMS'});
legend('Exact, frictionless', 'Exact friction');

% 5. Plots M/T for the case with friction
figure;
plot(alpha_a(1:end-5)*180/pi,1./T_over_M_friction_try1(1:end-5),'r-');
hold on;
plot(alpha_a(1:end-5)*180/pi,1./T_over_M_friction_try2(1:end-5),'k-');
hold on;
P5=gca;
xlabel('\alpha_a [deg]');
ylabel('M/T [-]');
title({'M/T, with friction: Ratio between net input and net output torques.';...
    'Exact (line) vs ADAMS envelope (dashed)'});

% 6. Plots F/M with and without friction
figure;
plot(alpha_a*180/pi,F_over_M,'g');hold on;
plot(alpha_a*180/pi,F_over_M_friction_try1,'r');hold on;
plot(alpha_a*180/pi,F_over_M_friction_try2,'k');
P6=gca;
xlabel('\alpha_a [deg]');
ylabel('F/M [1/mm]');
title({'F/M: Ratio between Body 2 force and net input torque.'});

% 7. Plots F/M with friction
figure;
plot(alpha_a(1:end-4)*180/pi,F_over_M_friction_try1(1:end-4),'k');hold on;
plot(alpha_a(1:end-4)*180/pi,F_over_M_friction_try2(1:end-4),'r');
P7=gca;
xlabel('\alpha_a [deg]');
ylabel('F/M [1/mm]');
title({'F/M, with friction: Ratio between Body 2 force and net input torque.';...
    'Exact (line) vs ADAMS envelope (dashed)'});

%--------------------------------------------------------------------------
% END
% PLOTTING THE EXACT RESULTS
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================
% NOTE:
% THE FOLLOWING LINES JUST PLOT THE ADAMS RESULTS FOR THE RSSR CASE SHOWN 
% IN  THE  PAPER.  IF YOU DONT WANT THESE PLOTS, DELETE ALL THE REMAINING 
% LINES.
%==========================================================================
%==========================================================================

% Declare global variables to recover ADAMS simulation results
global ADAMS_a ADAMS_c ADAMS_M ADAMS_T ADAMS_a_FRIC 
global ADAMS_T_FRIC ADAMS_M_FRIC ADAMS_dcsida
global ADAMS_a_FRIC_1 ADAMS_T_FRIC_1 ADAMS_M_FRIC_1
global ADAMS_a_FRIC_2 ADAMS_T_FRIC_2 ADAMS_M_FRIC_2
global ADAMS_a_FRIC_3 ADAMS_T_FRIC_3 ADAMS_M_FRIC_3
global ADAMS_a_FRIC_4 ADAMS_T_FRIC_4 ADAMS_M_FRIC_4
global ADAMS_F_OVER_M ADAMS_F_OVER_M_FRIC
global ADAMS_F_OVER_M_FRIC_1 ADAMS_F_OVER_M_FRIC_2
global ADAMS_F_OVER_M_FRIC_3 ADAMS_F_OVER_M_FRIC_4

%--------------------------------------------------------------------------
% START
% PLOTTING ADAMS RESULTS FOR THE APPLICATION OF THE PAPER 
%--------------------------------------------------------------------------

% Internal function callout to restore ADAMS simulation results.  They will
% be stored in the global variables defined in the beginning of the script.
restore_fcn;

% 1. Plots the input-output angles
axes(P1);
hold on;
plot(ADAMS_a(1:4:end)*180/pi,...
    ADAMS_c(1:4:end)*180/pi,'k.','markersize',10);
legend('Exact','ADAMS');

% 2. Plots normalized spin
axes(P2);
hold on;
plot(ADAMS_a_FRIC(1:4:end)*180/pi,...
    ADAMS_dcsida(1:4:end),'k.','markersize',10);
legend('Exact','ADAMS');

% 3. Plots pseudo-roll angle
    % 3.a Calculates the PSEUDO-roll angle:
    ADAMS_dcsida(1)=ADAMS_dcsida(2);
    csi_ADAMS=cumsum(ADAMS_dcsida.*gradient(ADAMS_a_FRIC));

    % 3.b. Applies a certain initial condition
    % (constant of integration - arbitrary)
    csi_ADAMS=csi_ADAMS-csi_ADAMS(1);
    csi_ADAMS=csi_ADAMS-CONST;
    
    % 3.c Plots pseudo-roll angle
    axes(P3);
    hold on;
    plot(ADAMS_a_FRIC(1:4:end)*180/pi,...
        csi_ADAMS(1:4:end)*180/pi,'k.','markersize',10);
    legend('Exact','ADAMS');

% 4. Plots M/T
axes(P4);
hold on;
plot(ADAMS_a(1:4:end)*180/pi,...
    ADAMS_M(1:4:end)/ADAMS_T,'g.','markersize',10);
hold on;
plot(ADAMS_a_FRIC(1:4:end-1)*180/pi,...
    ADAMS_M_FRIC(1:4:end-1)/ADAMS_T_FRIC,'r.','markersize',10);
legend('Exact, Frictionless','Exact, Friction','Exact, Friction',...
    'ADAMS, Friction', 'ADAMS, Frictionless');

% 5. Computes and plots the envelope of T/M with friction from ADAMS
%
% Notice that:
% - ADAMS friction models are different from those used here.
% - Thus, exact match between current exact and ADAMS sols is impossible.
%
% HOWEVER:
% - If all COF in ADAMS revolute joints were multipled by a COF correction 
%   factor 1/2<=lambda<=3/2, then the results must match.
% - Since this constant is variable and unknown beforehand, ENVELOPES were
%   computed from ADAMS.
% - Exact solution must lie within these ENVELOPES.
%
% - Since the RSSR has 2 revolute joints, 2 different lambdas were used.
%   The following cases were used to construct the envelopes:
%    lambda_i = 1/2, lambda_j = 1/2
%    lambda_i = 1/2, lambda_j = 3/2
%    lambda_i = 3/2, lambda_j = 1/2
%    lambda_i = 3/2, lambda_j = 3/2
%
% 5.1. First, all "y" must correspond to the same "x".
ADAMS_M_FRIC_1=interp1(ADAMS_a_FRIC_1,ADAMS_M_FRIC_1,ADAMS_a_FRIC_2);
ADAMS_M_FRIC_3=interp1(ADAMS_a_FRIC_3,ADAMS_M_FRIC_3,ADAMS_a_FRIC_2);
ADAMS_M_FRIC_4=interp1(ADAMS_a_FRIC_4,ADAMS_M_FRIC_4,ADAMS_a_FRIC_2);

% 5.2. Then, the envelope (max and min) is computed.
ADAMS_ENV_MAX=max([ADAMS_M_FRIC_1/ADAMS_T_FRIC_1,...
    ADAMS_M_FRIC_2/ADAMS_T_FRIC_2,...
    ADAMS_M_FRIC_3/ADAMS_T_FRIC_3,...
    ADAMS_M_FRIC_4/ADAMS_T_FRIC_4],[],2);
ADAMS_ENV_MIN=min([ADAMS_M_FRIC_1/ADAMS_T_FRIC_1,...
    ADAMS_M_FRIC_2/ADAMS_T_FRIC_2,...
    ADAMS_M_FRIC_3/ADAMS_T_FRIC_3,...
    ADAMS_M_FRIC_4/ADAMS_T_FRIC_4],[],2);

% 5.3. Finally, the envelope is plotted
axes(P5);
hold on;
plot(ADAMS_a_FRIC_2(1:4:end-4)*180/pi,...
    ADAMS_ENV_MAX(1:4:end-4),'k--','markersize',10);
hold on;
plot(ADAMS_a_FRIC_2(1:4:end-4)*180/pi,...
    ADAMS_ENV_MIN(1:4:end-4),'k--','markersize',10);

% 6. Plots F/M with and without friction
axes(P6);
hold on;
plot(ADAMS_a(1:6:end)*180/pi,...
    ADAMS_F_OVER_M(1:6:end),'k.','markersize',10);
hold on;

% Creates  a  filter  to control the density of points.
% Close to singularity they tend to be too spaced, and
% this  filter  increases the plot density there. This
% creates a more harmonic figure.
FILTER=ADAMS_a_FRIC*180/pi<-40 & ADAMS_a_FRIC*180/pi>-55;
FILTER=(FILTER.*1)./(FILTER.*1); 
plot(ADAMS_a_FRIC(1:3:end)*180/pi.*FILTER(1:3:end),...
    ADAMS_F_OVER_M_FRIC(1:3:end).*FILTER(1:3:end),'k.','markersize',10);
hold on;
plot(ADAMS_a_FRIC(3:6:end-4)*180/pi,...
    ADAMS_F_OVER_M_FRIC(3:6:end-4),'k.','markersize',10);
hold on;
plot(ADAMS_a_FRIC*180/pi.*FILTER,...
    ADAMS_F_OVER_M_FRIC.*FILTER,'k.','markersize',10);
set(gca,'ylim',[-.04,.04]);
legend('Exact, Frictionless', 'Exact, Friction','Exact, Friction','ADAMS');


% 7. Computes and plots the envelope of F/M with friction from ADAMS
%
% Notice that:
% - ADAMS friction models are different from those used here.
% - Thus, exact match between current exact and ADAMS sols is impossible.
%
% HOWEVER:
% - If all COF in ADAMS revolute joints were multipled by a COF correction 
%   factor 1/2<=lambda<=3/2, then the results must match.
% - Since this constant is variable and unknown beforehand, ENVELOPES were
%   computed from ADAMS.
% - Exact solution must lie within these ENVELOPES.
%
% - Since the RSSR has 2 revolute joints, 2 different lambdas were used.
%   The following cases were used to construct the envelopes:
%    lambda_i = 1/2, lambda_j = 1/2
%    lambda_i = 1/2, lambda_j = 3/2
%    lambda_i = 3/2, lambda_j = 1/2
%    lambda_i = 3/2, lambda_j = 3/2
%
% 7.1. First, all "y" must correspond to the same "x".
ADAMS_F_OVER_M_FRIC_1=interp1(ADAMS_a_FRIC_1,ADAMS_F_OVER_M_FRIC_1,ADAMS_a_FRIC_2);
ADAMS_F_OVER_M_FRIC_3=interp1(ADAMS_a_FRIC_3,ADAMS_F_OVER_M_FRIC_3,ADAMS_a_FRIC_2);
ADAMS_F_OVER_M_FRIC_4=interp1(ADAMS_a_FRIC_4,ADAMS_F_OVER_M_FRIC_4,ADAMS_a_FRIC_2);

% 7.2. Then, the envelope (max and min) is computed.
ADAMS_ENV_MAX=max([ADAMS_F_OVER_M_FRIC_1,...
    ADAMS_F_OVER_M_FRIC_2,...
    ADAMS_F_OVER_M_FRIC_3,...
    ADAMS_F_OVER_M_FRIC_4],[],2);
ADAMS_ENV_MIN=min([ADAMS_F_OVER_M_FRIC_1,...
    ADAMS_F_OVER_M_FRIC_2,...
    ADAMS_F_OVER_M_FRIC_3,...
    ADAMS_F_OVER_M_FRIC_4],[],2);

% Creates  a  filter close to singularity to avoid noisy lines
FILTER= ADAMS_a_FRIC_2*180/pi>-42.97 | ADAMS_a_FRIC_2*180/pi<-52.02;
FILTER=(FILTER.*1)./(FILTER.*1);
axes(P7);
hold on;
plot(ADAMS_a_FRIC_2.*FILTER.*180/pi,...
    ADAMS_ENV_MAX.*FILTER,'k--','markersize',10);
hold on;
plot(ADAMS_a_FRIC_2.*FILTER.*180/pi,...
    ADAMS_ENV_MIN.*FILTER,'k--','markersize',10);
set(gca,'ylim',[-.04,.04]);
%--------------------------------------------------------------------------
% END
% PLOTTING ADAMS RESULTS FOR THE APPLICATION OF THE PAPER
%--------------------------------------------------------------------------


% IF  YOU  DONT  WANT  TO USE THE ADAMS RESULTS SHOWN IN THE PAPER, YOU CAN 
% DELETE THE FOLLOWING FUNCTION ENTIRELY.

function restore_fcn

% This internal function just populate the following global variables with 
% results from ADAMS simulation.
% 
% The article compared exact solution with ADAMS. The RSSR mechanism of the 
% paper was simulated in ADAMS and the results were stored in the sequence.
%
% When this function is called out, those ADAMS results are stored in the
% global variables, which can then be used to generate the graphs shown in
% the paper.

global ADAMS_a ADAMS_c ADAMS_M ADAMS_T ADAMS_a_FRIC 
global ADAMS_T_FRIC ADAMS_M_FRIC ADAMS_dcsida
global ADAMS_a_FRIC_1 ADAMS_T_FRIC_1 ADAMS_M_FRIC_1
global ADAMS_a_FRIC_2 ADAMS_T_FRIC_2 ADAMS_M_FRIC_2
global ADAMS_a_FRIC_3 ADAMS_T_FRIC_3 ADAMS_M_FRIC_3
global ADAMS_a_FRIC_4 ADAMS_T_FRIC_4 ADAMS_M_FRIC_4 
global ADAMS_F_OVER_M ADAMS_F_OVER_M_FRIC
global ADAMS_F_OVER_M_FRIC_1 ADAMS_F_OVER_M_FRIC_2
global ADAMS_F_OVER_M_FRIC_3 ADAMS_F_OVER_M_FRIC_4

ADAMS_a=[         0
   -0.0175
   -0.0349
   -0.0524
   -0.0698
   -0.0873
   -0.1050
   -0.1220
   -0.1400
   -0.1570
   -0.1750
   -0.1920
   -0.2090
   -0.2270
   -0.2440
   -0.2620
   -0.2790
   -0.2970
   -0.3140
   -0.3320
   -0.3490
   -0.3670
   -0.3840
   -0.4010
   -0.4190
   -0.4360
   -0.4540
   -0.4710
   -0.4890
   -0.5060
   -0.5240
   -0.5410
   -0.5590
   -0.5760
   -0.5930
   -0.6110
   -0.6280
   -0.6460
   -0.6630
   -0.6810
   -0.6980
   -0.7160
   -0.7330
   -0.7500
   -0.7680
   -0.7850
   -0.8030
   -0.8200
   -0.8380
   -0.8550
   -0.8730
   -0.8900
   -0.9080
   -0.9250
   -0.9420
   -0.9600
   -0.9770
   -0.9950
   -1.0100
   -1.0300
   -1.0500
   -1.0600
   -1.0800
   -1.1000
   -1.1200
   -1.1300
   -1.1500
   -1.1700
   -1.1900
   -1.2000
   -1.2200
   -1.2400
   -1.2600
   -1.2700
   -1.2900
   -1.3100
   -1.3300
   -1.3400
   -1.3600
   -1.3800
   -1.4000
   -1.4100
   -1.4300
   -1.4500
   -1.4700
   -1.4800
   -1.5000
   -1.5200
   -1.5400
   -1.5500
   -1.5700
   -1.5900
   -1.6100
   -1.6200
   -1.6400
   -1.6600
   -1.6800
   -1.6900
   -1.7100
   -1.7300
   -1.7500
   -1.7600
   -1.7800
   -1.8000
   -1.8200
   -1.8300
   -1.8500
   -1.8700
   -1.8800
   -1.9000
   -1.9200
   -1.9400
   -1.9500
   -1.9700
   -1.9900
   -2.0100
   -2.0200
   -2.0400
   -2.0600
   -2.0800
   -2.0900
   -2.1100
   -2.1300
   -2.1500
   -2.1600
   -2.1800
   -2.2000
   -2.2200
   -2.2300
   -2.2500
   -2.2700
   -2.2900
   -2.3000
   -2.3200
   -2.3400
   -2.3600
   -2.3700
   -2.3900
   -2.4100
   -2.4300
   -2.4400
   -2.4600
   -2.4800
   -2.5000
   -2.5100
   -2.5300
   -2.5500
   -2.5700
   -2.5800
   -2.6000
   -2.6200
   -2.6400
   -2.6500
   -2.6700
   -2.6900
   -2.7100
   -2.7200
   -2.7400
   -2.7600
   -2.7800
   -2.7900
   -2.8100
   -2.8300
   -2.8400
   -2.8600
   -2.8800
   -2.9000
   -2.9100
   -2.9300
   -2.9500
   -2.9700
   -2.9800
   -3.0000
   -3.0200
   -3.0400
   -3.0500
   -3.0700
   -3.0900
   -3.1100
   -3.1200
   -3.1400
   -3.1600
   -3.1800
   -3.1900
   -3.2100
   -3.2300
   -3.2500
   -3.2600
   -3.2800
   -3.3000
   -3.3200
   -3.3300
   -3.3500
   -3.3700
   -3.3900
   -3.4000
   -3.4200
   -3.4400
   -3.4600
   -3.4700
   -3.4900
   -3.5100
   -3.5300
   -3.5400
   -3.5600
   -3.5800
   -3.6000
   -3.6100
   -3.6300
   -3.6500
   -3.6700
   -3.6800
   -3.7000
   -3.7200
   -3.7400
   -3.7500
   -3.7700
   -3.7900
   -3.8000
   -3.8200
   -3.8400
   -3.8600
   -3.8700
   -3.8900
   -3.9100
   -3.9300
   -3.9400
   -3.9600
   -3.9800
   -4.0000
   -4.0100
   -4.0300
   -4.0500
   -4.0700
   -4.0800
   -4.1000
   -4.1200
   -4.1400];

ADAMS_c=[    1.7100
    1.6900
    1.6700
    1.6500
    1.6300
    1.6100
    1.5900
    1.5800
    1.5600
    1.5400
    1.5300
    1.5100
    1.5000
    1.4800
    1.4700
    1.4600
    1.4400
    1.4300
    1.4200
    1.4000
    1.3900
    1.3800
    1.3700
    1.3600
    1.3500
    1.3400
    1.3300
    1.3200
    1.3100
    1.3000
    1.2900
    1.2800
    1.2700
    1.2600
    1.2600
    1.2500
    1.2400
    1.2300
    1.2300
    1.2200
    1.2200
    1.2100
    1.2100
    1.2000
    1.2000
    1.1900
    1.1900
    1.1800
    1.1800
    1.1700
    1.1700
    1.1700
    1.1700
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1600
    1.1700
    1.1700
    1.1700
    1.1700
    1.1800
    1.1800
    1.1900
    1.1900
    1.1900
    1.2000
    1.2000
    1.2100
    1.2200
    1.2200
    1.2300
    1.2400
    1.2400
    1.2500
    1.2600
    1.2600
    1.2700
    1.2800
    1.2900
    1.3000
    1.3100
    1.3200
    1.3300
    1.3400
    1.3500
    1.3600
    1.3700
    1.3800
    1.3900
    1.4000
    1.4100
    1.4200
    1.4300
    1.4500
    1.4600
    1.4700
    1.4800
    1.5000
    1.5100
    1.5200
    1.5400
    1.5500
    1.5600
    1.5700
    1.5900
    1.6000
    1.6200
    1.6300
    1.6400
    1.6600
    1.6700
    1.6900
    1.7000
    1.7100
    1.7300
    1.7400
    1.7600
    1.7700
    1.7800
    1.8000
    1.8100
    1.8300
    1.8400
    1.8600
    1.8700
    1.8800
    1.9000
    1.9100
    1.9300
    1.9400
    1.9500
    1.9700
    1.9800
    1.9900
    2.0100
    2.0200
    2.0300
    2.0500
    2.0600
    2.0700
    2.0900
    2.1000
    2.1100
    2.1300
    2.1400
    2.1500
    2.1700
    2.1800
    2.1900
    2.2000
    2.2200
    2.2300
    2.2400
    2.2500
    2.2700
    2.2800
    2.2900
    2.3000
    2.3100
    2.3300
    2.3400
    2.3500
    2.3600
    2.3700
    2.3900
    2.4000
    2.4100
    2.4200
    2.4300
    2.4400
    2.4500
    2.4700
    2.4800
    2.4900
    2.5000
    2.5100
    2.5200
    2.5300
    2.5400
    2.5600
    2.5700
    2.5800
    2.5900
    2.6000
    2.6100
    2.6200
    2.6300
    2.6500
    2.6600
    2.6700
    2.6800
    2.6900
    2.7000
    2.7100
    2.7300
    2.7400
    2.7500
    2.7600
    2.7700
    2.7900
    2.8000
    2.8100
    2.8200
    2.8400
    2.8500
    2.8600
    2.8800
    2.8900
    2.9000
    2.9200
    2.9300
    2.9500
    2.9600
    2.9800
    3.0000
    3.0200
    3.0300
    3.0500
    3.0800
    3.1000
    3.1200
    3.1500
    3.1800
    3.2200
    3.2700];

ADAMS_T=-1;

ADAMS_M=[    1.1600
    1.1300
    1.1000
    1.0700
    1.0400
    1.0200
    0.9910
    0.9650
    0.9400
    0.9160
    0.8930
    0.8700
    0.8470
    0.8250
    0.8040
    0.7830
    0.7620
    0.7420
    0.7220
    0.7020
    0.6830
    0.6640
    0.6450
    0.6260
    0.6080
    0.5900
    0.5720
    0.5540
    0.5370
    0.5190
    0.5020
    0.4850
    0.4680
    0.4510
    0.4340
    0.4170
    0.4010
    0.3840
    0.3680
    0.3510
    0.3350
    0.3190
    0.3030
    0.2860
    0.2700
    0.2540
    0.2380
    0.2220
    0.2050
    0.1890
    0.1730
    0.1570
    0.1410
    0.1250
    0.1080
    0.0920
    0.0758
    0.0594
    0.0431
    0.0267
    0.0103
   -0.0062
   -0.0227
   -0.0392
   -0.0558
   -0.0724
   -0.0891
   -0.1060
   -0.1220
   -0.1390
   -0.1560
   -0.1730
   -0.1890
   -0.2060
   -0.2230
   -0.2400
   -0.2560
   -0.2730
   -0.2900
   -0.3060
   -0.3230
   -0.3400
   -0.3560
   -0.3720
   -0.3890
   -0.4050
   -0.4210
   -0.4370
   -0.4520
   -0.4680
   -0.4830
   -0.4980
   -0.5130
   -0.5280
   -0.5420
   -0.5570
   -0.5700
   -0.5840
   -0.5970
   -0.6100
   -0.6230
   -0.6350
   -0.6470
   -0.6590
   -0.6700
   -0.6810
   -0.6910
   -0.7010
   -0.7110
   -0.7200
   -0.7290
   -0.7370
   -0.7450
   -0.7520
   -0.7590
   -0.7660
   -0.7720
   -0.7770
   -0.7820
   -0.7870
   -0.7910
   -0.7950
   -0.7990
   -0.8020
   -0.8040
   -0.8070
   -0.8080
   -0.8100
   -0.8110
   -0.8120
   -0.8120
   -0.8120
   -0.8120
   -0.8110
   -0.8110
   -0.8090
   -0.8080
   -0.8060
   -0.8050
   -0.8020
   -0.8000
   -0.7980
   -0.7950
   -0.7920
   -0.7890
   -0.7860
   -0.7820
   -0.7790
   -0.7760
   -0.7720
   -0.7680
   -0.7640
   -0.7600
   -0.7570
   -0.7530
   -0.7490
   -0.7450
   -0.7400
   -0.7360
   -0.7320
   -0.7280
   -0.7240
   -0.7200
   -0.7160
   -0.7120
   -0.7080
   -0.7050
   -0.7010
   -0.6970
   -0.6930
   -0.6900
   -0.6860
   -0.6830
   -0.6790
   -0.6760
   -0.6730
   -0.6700
   -0.6670
   -0.6640
   -0.6620
   -0.6590
   -0.6570
   -0.6550
   -0.6530
   -0.6510
   -0.6490
   -0.6470
   -0.6460
   -0.6440
   -0.6430
   -0.6420
   -0.6420
   -0.6410
   -0.6410
   -0.6410
   -0.6410
   -0.6420
   -0.6420
   -0.6440
   -0.6450
   -0.6460
   -0.6480
   -0.6510
   -0.6540
   -0.6570
   -0.6600
   -0.6640
   -0.6690
   -0.6740
   -0.6790
   -0.6850
   -0.6920
   -0.7000
   -0.7080
   -0.7180
   -0.7280
   -0.7390
   -0.7520
   -0.7660
   -0.7820
   -0.7990
   -0.8180
   -0.8400
   -0.8640
   -0.8920
   -0.9230
   -0.9600
   -1.0000
   -1.0500
   -1.1100
   -1.1800
   -1.2600
   -1.3700
   -1.5100
   -1.7000
   -1.9900
   -2.4800
   -3.5900];

ADAMS_F_OVER_M=[6.77E-03
6.82E-03
6.87E-03
6.92E-03
6.98E-03
7.05E-03
7.12E-03
7.20E-03
7.28E-03
7.37E-03
7.47E-03
7.57E-03
7.67E-03
7.79E-03
7.91E-03
8.04E-03
8.18E-03
8.32E-03
8.47E-03
8.64E-03
8.81E-03
8.99E-03
9.19E-03
9.39E-03
9.61E-03
9.84E-03
1.01E-02
1.04E-02
1.06E-02
1.09E-02
1.13E-02
1.16E-02
1.20E-02
1.24E-02
1.28E-02
1.33E-02
1.38E-02
1.43E-02
1.49E-02
1.56E-02
1.63E-02
1.71E-02
1.80E-02
1.90E-02
2.01E-02
2.13E-02
2.27E-02
2.44E-02
2.63E-02
2.85E-02
3.11E-02
3.43E-02
3.83E-02
4.32E-02
4.97E-02
5.85E-02
7.11E-02
9.06E-02
1.25E-01
2.02E-01
5.26E-01
-8.70E-01
-2.38E-01
-1.38E-01
-9.70E-02
-7.48E-02
-6.09E-02
-5.14E-02
-4.44E-02
-3.91E-02
-3.50E-02
-3.16E-02
-2.89E-02
-2.66E-02
-2.46E-02
-2.29E-02
-2.14E-02
-2.02E-02
-1.90E-02
-1.80E-02
-1.71E-02
-1.63E-02
-1.56E-02
-1.49E-02
-1.43E-02
-1.38E-02
-1.33E-02
-1.28E-02
-1.24E-02
-1.20E-02
-1.16E-02
-1.12E-02
-1.09E-02
-1.06E-02
-1.04E-02
-1.01E-02
-9.86E-03
-9.63E-03
-9.42E-03
-9.22E-03
-9.03E-03
-8.86E-03
-8.70E-03
-8.54E-03
-8.40E-03
-8.26E-03
-8.13E-03
-8.01E-03
-7.90E-03
-7.79E-03
-7.69E-03
-7.60E-03
-7.51E-03
-7.43E-03
-7.35E-03
-7.28E-03
-7.21E-03
-7.15E-03
-7.09E-03
-7.04E-03
-6.99E-03
-6.94E-03
-6.90E-03
-6.86E-03
-6.82E-03
-6.79E-03
-6.76E-03
-6.73E-03
-6.71E-03
-6.69E-03
-6.67E-03
-6.66E-03
-6.64E-03
-6.63E-03
-6.62E-03
-6.62E-03
-6.62E-03
-6.61E-03
-6.61E-03
-6.62E-03
-6.62E-03
-6.63E-03
-6.64E-03
-6.65E-03
-6.66E-03
-6.67E-03
-6.69E-03
-6.71E-03
-6.72E-03
-6.74E-03
-6.77E-03
-6.79E-03
-6.81E-03
-6.84E-03
-6.87E-03
-6.90E-03
-6.93E-03
-6.96E-03
-6.99E-03
-7.03E-03
-7.06E-03
-7.10E-03
-7.14E-03
-7.17E-03
-7.21E-03
-7.26E-03
-7.30E-03
-7.34E-03
-7.38E-03
-7.43E-03
-7.48E-03
-7.52E-03
-7.57E-03
-7.62E-03
-7.67E-03
-7.72E-03
-7.77E-03
-7.82E-03
-7.88E-03
-7.93E-03
-7.99E-03
-8.04E-03
-8.10E-03
-8.16E-03
-8.21E-03
-8.27E-03
-8.33E-03
-8.39E-03
-8.45E-03
-8.51E-03
-8.57E-03
-8.63E-03
-8.69E-03
-8.76E-03
-8.82E-03
-8.88E-03
-8.95E-03
-9.01E-03
-9.07E-03
-9.14E-03
-9.20E-03
-9.26E-03
-9.33E-03
-9.39E-03
-9.45E-03
-9.51E-03
-9.58E-03
-9.64E-03
-9.70E-03
-9.76E-03
-9.82E-03
-9.88E-03
-9.93E-03
-9.99E-03
-1.00E-02
-1.01E-02
-1.01E-02
-1.02E-02
-1.02E-02
-1.03E-02
-1.03E-02
-1.04E-02
-1.04E-02
-1.04E-02
-1.04E-02
-1.05E-02
-1.05E-02
-1.05E-02
-1.05E-02
-1.05E-02
-1.04E-02
-1.04E-02
-1.04E-02
-1.03E-02
-1.02E-02
-1.01E-02
-9.97E-03
-9.74E-03
];

ADAMS_a_FRIC=[0.00E+00
-1.75E-02
-3.49E-02
-5.24E-02
-6.98E-02
-8.73E-02
-1.05E-01
-1.22E-01
-1.40E-01
-1.57E-01
-1.75E-01
-1.92E-01
-2.09E-01
-2.27E-01
-2.44E-01
-2.62E-01
-2.79E-01
-2.97E-01
-3.14E-01
-3.32E-01
-3.49E-01
-3.67E-01
-3.84E-01
-4.01E-01
-4.19E-01
-4.36E-01
-4.54E-01
-4.71E-01
-4.89E-01
-5.06E-01
-5.24E-01
-5.41E-01
-5.59E-01
-5.76E-01
-5.93E-01
-6.11E-01
-6.28E-01
-6.46E-01
-6.63E-01
-6.81E-01
-6.98E-01
-7.16E-01
-7.33E-01
-7.50E-01
-7.68E-01
-7.85E-01
-8.03E-01
-8.20E-01
-8.38E-01
-8.55E-01
-8.73E-01
-8.90E-01
-9.08E-01
-9.25E-01
-9.42E-01
-9.60E-01
-9.77E-01
-9.95E-01
-1.01E+00
-1.03E+00
-1.05E+00
-1.06E+00
-1.08E+00
-1.10E+00
-1.12E+00
-1.13E+00
-1.15E+00
-1.17E+00
-1.19E+00
-1.20E+00
-1.22E+00
-1.24E+00
-1.26E+00
-1.27E+00
-1.29E+00
-1.31E+00
-1.33E+00
-1.34E+00
-1.36E+00
-1.38E+00
-1.40E+00
-1.41E+00
-1.43E+00
-1.45E+00
-1.47E+00
-1.48E+00
-1.50E+00
-1.52E+00
-1.54E+00
-1.55E+00
-1.57E+00
-1.59E+00
-1.61E+00
-1.62E+00
-1.64E+00
-1.66E+00
-1.68E+00
-1.69E+00
-1.71E+00
-1.73E+00
-1.75E+00
-1.76E+00
-1.78E+00
-1.80E+00
-1.82E+00
-1.83E+00
-1.85E+00
-1.87E+00
-1.88E+00
-1.90E+00
-1.92E+00
-1.94E+00
-1.95E+00
-1.97E+00
-1.99E+00
-2.01E+00
-2.02E+00
-2.04E+00
-2.06E+00
-2.08E+00
-2.09E+00
-2.11E+00
-2.13E+00
-2.15E+00
-2.16E+00
-2.18E+00
-2.20E+00
-2.22E+00
-2.23E+00
-2.25E+00
-2.27E+00
-2.29E+00
-2.30E+00
-2.32E+00
-2.34E+00
-2.36E+00
-2.37E+00
-2.39E+00
-2.41E+00
-2.43E+00
-2.44E+00
-2.46E+00
-2.48E+00
-2.50E+00
-2.51E+00
-2.53E+00
-2.55E+00
-2.57E+00
-2.58E+00
-2.60E+00
-2.62E+00
-2.64E+00
-2.65E+00
-2.67E+00
-2.69E+00
-2.71E+00
-2.72E+00
-2.74E+00
-2.76E+00
-2.78E+00
-2.79E+00
-2.81E+00
-2.83E+00
-2.84E+00
-2.86E+00
-2.88E+00
-2.90E+00
-2.91E+00
-2.93E+00
-2.95E+00
-2.97E+00
-2.98E+00
-3.00E+00
-3.02E+00
-3.04E+00
-3.05E+00
-3.07E+00
-3.09E+00
-3.11E+00
-3.12E+00
-3.14E+00
-3.16E+00
-3.18E+00
-3.19E+00
-3.21E+00
-3.23E+00
-3.25E+00
-3.26E+00
-3.28E+00
-3.30E+00
-3.32E+00
-3.33E+00
-3.35E+00
-3.37E+00
-3.39E+00
-3.40E+00
-3.42E+00
-3.44E+00
-3.46E+00
-3.47E+00
-3.49E+00
-3.51E+00
-3.53E+00
-3.54E+00
-3.56E+00
-3.58E+00
-3.60E+00
-3.61E+00
-3.63E+00
-3.65E+00
-3.67E+00
-3.68E+00
-3.70E+00
-3.72E+00
-3.74E+00
-3.75E+00
-3.77E+00
-3.79E+00
-3.80E+00
-3.82E+00
-3.84E+00
-3.86E+00
-3.87E+00
-3.89E+00
-3.91E+00
-3.93E+00
-3.94E+00
-3.96E+00
-3.98E+00
-4.00E+00
-4.01E+00
-4.03E+00
-4.05E+00
-4.07E+00
-4.08E+00
-4.10E+00
];

ADAMS_T_FRIC=-50;

ADAMS_M_FRIC=[3.70E+01
3.63E+01
3.53E+01
3.44E+01
3.35E+01
3.25E+01
3.16E+01
3.08E+01
2.99E+01
2.90E+01
2.82E+01
2.74E+01
2.65E+01
2.57E+01
2.49E+01
2.41E+01
2.33E+01
2.25E+01
2.17E+01
2.10E+01
2.02E+01
1.94E+01
1.87E+01
1.79E+01
1.72E+01
1.64E+01
1.57E+01
1.50E+01
1.42E+01
1.35E+01
1.28E+01
1.20E+01
1.13E+01
1.06E+01
9.85E+00
9.13E+00
8.40E+00
7.68E+00
6.96E+00
6.23E+00
5.51E+00
4.78E+00
4.05E+00
3.33E+00
2.60E+00
1.87E+00
1.13E+00
4.01E-01
-3.35E-01
-1.07E+00
-1.81E+00
-2.56E+00
-3.30E+00
-4.05E+00
-4.80E+00
-5.55E+00
-6.31E+00
-7.07E+00
-7.83E+00
-8.59E+00
-9.36E+00
-1.11E+01
-1.19E+01
-1.27E+01
-1.36E+01
-1.44E+01
-1.52E+01
-1.61E+01
-1.69E+01
-1.77E+01
-1.86E+01
-1.94E+01
-2.02E+01
-2.11E+01
-2.19E+01
-2.27E+01
-2.35E+01
-2.44E+01
-2.52E+01
-2.60E+01
-2.68E+01
-2.76E+01
-2.84E+01
-2.92E+01
-3.00E+01
-3.08E+01
-3.16E+01
-3.24E+01
-3.31E+01
-3.39E+01
-3.46E+01
-3.54E+01
-3.61E+01
-3.68E+01
-3.75E+01
-3.82E+01
-3.89E+01
-3.96E+01
-4.02E+01
-4.09E+01
-4.15E+01
-4.22E+01
-4.28E+01
-4.35E+01
-4.41E+01
-4.47E+01
-4.52E+01
-4.58E+01
-4.62E+01
-4.67E+01
-4.71E+01
-4.75E+01
-4.79E+01
-4.82E+01
-4.85E+01
-4.88E+01
-4.91E+01
-4.93E+01
-4.96E+01
-4.98E+01
-5.00E+01
-5.01E+01
-5.03E+01
-5.04E+01
-5.05E+01
-5.06E+01
-5.07E+01
-5.07E+01
-5.08E+01
-5.08E+01
-5.08E+01
-5.07E+01
-5.07E+01
-5.07E+01
-5.06E+01
-5.05E+01
-5.04E+01
-5.03E+01
-5.02E+01
-5.01E+01
-5.00E+01
-4.98E+01
-4.97E+01
-4.95E+01
-4.94E+01
-4.92E+01
-4.90E+01
-4.89E+01
-4.87E+01
-4.85E+01
-4.83E+01
-4.81E+01
-4.79E+01
-4.78E+01
-4.76E+01
-4.74E+01
-4.72E+01
-4.70E+01
-4.68E+01
-4.66E+01
-4.64E+01
-4.62E+01
-4.61E+01
-4.59E+01
-4.57E+01
-4.56E+01
-4.54E+01
-4.53E+01
-4.51E+01
-4.50E+01
-4.48E+01
-4.47E+01
-4.46E+01
-4.45E+01
-4.44E+01
-4.43E+01
-4.42E+01
-4.41E+01
-4.41E+01
-4.40E+01
-4.40E+01
-4.39E+01
-4.39E+01
-4.39E+01
-4.39E+01
-4.39E+01
-4.40E+01
-4.40E+01
-4.41E+01
-4.42E+01
-4.43E+01
-4.44E+01
-4.45E+01
-4.47E+01
-4.49E+01
-4.51E+01
-4.53E+01
-4.56E+01
-4.58E+01
-4.62E+01
-4.65E+01
-4.69E+01
-4.73E+01
-4.78E+01
-4.83E+01
-4.88E+01
-4.94E+01
-5.00E+01
-5.08E+01
-5.15E+01
-5.24E+01
-5.33E+01
-5.43E+01
-5.55E+01
-5.67E+01
-5.81E+01
-5.96E+01
-6.13E+01
-6.32E+01
-6.53E+01
-6.77E+01
-7.05E+01
-7.36E+01
-7.73E+01
-8.17E+01
-8.68E+01
-9.30E+01
-1.00E+02
-1.10E+02
-1.21E+02
-1.37E+02
-1.58E+02
-1.90E+02
-2.42E+02
-3.45E+02
-6.58E+02
];

ADAMS_dcsida=[9.39E-13
6.36E-01
6.07E-01
5.78E-01
5.51E-01
5.24E-01
4.99E-01
4.74E-01
4.50E-01
4.27E-01
4.05E-01
3.84E-01
3.63E-01
3.43E-01
3.23E-01
3.04E-01
2.86E-01
2.68E-01
2.51E-01
2.34E-01
2.18E-01
2.02E-01
1.87E-01
1.72E-01
1.58E-01
1.44E-01
1.31E-01
1.18E-01
1.06E-01
9.35E-02
8.19E-02
7.07E-02
5.99E-02
4.94E-02
3.94E-02
2.98E-02
2.05E-02
1.16E-02
3.10E-03
-5.06E-03
-1.29E-02
-2.03E-02
-2.74E-02
-3.42E-02
-4.06E-02
-4.67E-02
-5.25E-02
-5.79E-02
-6.30E-02
-6.78E-02
-7.23E-02
-7.65E-02
-8.03E-02
-8.39E-02
-8.71E-02
-9.01E-02
-9.27E-02
-9.51E-02
-9.72E-02
-9.90E-02
-1.00E-01
-1.02E-01
-1.03E-01
-1.03E-01
-1.04E-01
-1.04E-01
-1.04E-01
-1.04E-01
-1.03E-01
-1.02E-01
-1.01E-01
-1.00E-01
-9.85E-02
-9.68E-02
-9.49E-02
-9.28E-02
-9.05E-02
-8.80E-02
-8.54E-02
-8.26E-02
-7.96E-02
-7.65E-02
-7.32E-02
-6.98E-02
-6.63E-02
-6.26E-02
-5.89E-02
-5.51E-02
-5.12E-02
-4.73E-02
-4.33E-02
-3.92E-02
-3.52E-02
-3.11E-02
-2.71E-02
-2.31E-02
-1.91E-02
-1.52E-02
-1.14E-02
-7.77E-03
-4.29E-03
-1.02E-03
1.98E-03
4.65E-03
6.93E-03
8.77E-03
1.01E-02
1.10E-02
1.13E-02
1.12E-02
1.08E-02
9.97E-03
8.88E-03
7.53E-03
5.97E-03
4.20E-03
2.24E-03
1.11E-04
-2.18E-03
-4.63E-03
-7.23E-03
-9.96E-03
-1.28E-02
-1.58E-02
-1.90E-02
-2.22E-02
-2.56E-02
-2.91E-02
-3.27E-02
-3.64E-02
-4.01E-02
-4.40E-02
-4.80E-02
-5.20E-02
-5.62E-02
-6.04E-02
-6.46E-02
-6.89E-02
-7.33E-02
-7.78E-02
-8.22E-02
-8.67E-02
-9.13E-02
-9.59E-02
-1.00E-01
-1.05E-01
-1.10E-01
-1.14E-01
-1.19E-01
-1.24E-01
-1.28E-01
-1.33E-01
-1.38E-01
-1.42E-01
-1.47E-01
-1.52E-01
-1.56E-01
-1.61E-01
-1.65E-01
-1.70E-01
-1.74E-01
-1.79E-01
-1.83E-01
-1.88E-01
-1.92E-01
-1.96E-01
-2.01E-01
-2.05E-01
-2.09E-01
-2.13E-01
-2.18E-01
-2.22E-01
-2.26E-01
-2.30E-01
-2.34E-01
-2.38E-01
-2.42E-01
-2.46E-01
-2.50E-01
-2.54E-01
-2.57E-01
-2.61E-01
-2.65E-01
-2.69E-01
-2.72E-01
-2.76E-01
-2.80E-01
-2.83E-01
-2.87E-01
-2.91E-01
-2.94E-01
-2.98E-01
-3.02E-01
-3.05E-01
-3.09E-01
-3.12E-01
-3.16E-01
-3.20E-01
-3.23E-01
-3.27E-01
-3.30E-01
-3.34E-01
-3.38E-01
-3.41E-01
-3.45E-01
-3.49E-01
-3.52E-01
-3.56E-01
-3.60E-01
-3.64E-01
-3.67E-01
-3.71E-01
-3.75E-01
-3.78E-01
-3.82E-01
-3.85E-01
-3.89E-01
-3.92E-01
-3.95E-01
-3.98E-01
-4.01E-01
-4.05E-01
-4.08E-01
-4.13E-01
-4.20E-01
-4.32E-01
-4.50E-01
-4.77E-01
-5.12E-01
-5.55E-01
-6.06E-01
-6.69E-01
-7.48E-01
-8.49E-01
-9.85E-01
-1.19E+00
];

ADAMS_a_FRIC_1=[0.00E+00
-1.75E-02
-3.49E-02
-5.24E-02
-6.98E-02
-8.73E-02
-1.05E-01
-1.22E-01
-1.40E-01
-1.57E-01
-1.75E-01
-1.92E-01
-2.09E-01
-2.27E-01
-2.44E-01
-2.62E-01
-2.79E-01
-2.97E-01
-3.14E-01
-3.32E-01
-3.49E-01
-3.67E-01
-3.84E-01
-4.01E-01
-4.19E-01
-4.36E-01
-4.54E-01
-4.71E-01
-4.89E-01
-5.06E-01
-5.24E-01
-5.41E-01
-5.59E-01
-5.76E-01
-5.93E-01
-6.11E-01
-6.28E-01
-6.46E-01
-6.63E-01
-6.81E-01
-6.98E-01
-7.16E-01
-7.33E-01
-7.50E-01
-7.68E-01
-7.85E-01
-8.03E-01
-8.20E-01
-8.38E-01
-8.55E-01
-8.73E-01
-8.90E-01
-9.08E-01
-9.25E-01
-9.42E-01
-9.60E-01
-9.77E-01
-9.95E-01
-1.01E+00
-1.03E+00
-1.05E+00
-1.06E+00
-1.08E+00
-1.10E+00
-1.12E+00
-1.13E+00
-1.15E+00
-1.17E+00
-1.19E+00
-1.20E+00
-1.22E+00
-1.24E+00
-1.26E+00
-1.27E+00
-1.29E+00
-1.31E+00
-1.33E+00
-1.34E+00
-1.36E+00
-1.38E+00
-1.40E+00
-1.41E+00
-1.43E+00
-1.45E+00
-1.47E+00
-1.48E+00
-1.50E+00
-1.52E+00
-1.54E+00
-1.55E+00
-1.57E+00
-1.59E+00
-1.61E+00
-1.62E+00
-1.64E+00
-1.66E+00
-1.68E+00
-1.69E+00
-1.71E+00
-1.73E+00
-1.75E+00
-1.76E+00
-1.78E+00
-1.80E+00
-1.82E+00
-1.83E+00
-1.85E+00
-1.87E+00
-1.88E+00
-1.90E+00
-1.92E+00
-1.94E+00
-1.95E+00
-1.97E+00
-1.99E+00
-2.01E+00
-2.02E+00
-2.04E+00
-2.06E+00
-2.08E+00
-2.09E+00
-2.11E+00
-2.13E+00
-2.15E+00
-2.16E+00
-2.18E+00
-2.20E+00
-2.22E+00
-2.23E+00
-2.25E+00
-2.27E+00
-2.29E+00
-2.30E+00
-2.32E+00
-2.34E+00
-2.36E+00
-2.37E+00
-2.39E+00
-2.41E+00
-2.43E+00
-2.44E+00
-2.46E+00
-2.48E+00
-2.50E+00
-2.51E+00
-2.53E+00
-2.55E+00
-2.57E+00
-2.58E+00
-2.60E+00
-2.62E+00
-2.64E+00
-2.65E+00
-2.67E+00
-2.69E+00
-2.71E+00
-2.72E+00
-2.74E+00
-2.76E+00
-2.78E+00
-2.79E+00
-2.81E+00
-2.83E+00
-2.84E+00
-2.86E+00
-2.88E+00
-2.90E+00
-2.91E+00
-2.93E+00
-2.95E+00
-2.97E+00
-2.98E+00
-3.00E+00
-3.02E+00
-3.04E+00
-3.05E+00
-3.07E+00
-3.09E+00
-3.11E+00
-3.12E+00
-3.14E+00
-3.16E+00
-3.18E+00
-3.19E+00
-3.21E+00
-3.23E+00
-3.25E+00
-3.26E+00
-3.28E+00
-3.30E+00
-3.32E+00
-3.33E+00
-3.35E+00
-3.37E+00
-3.39E+00
-3.40E+00
-3.42E+00
-3.44E+00
-3.46E+00
-3.47E+00
-3.49E+00
-3.51E+00
-3.53E+00
-3.54E+00
-3.56E+00
-3.58E+00
-3.60E+00
-3.61E+00
-3.63E+00
-3.65E+00
-3.67E+00
-3.68E+00
-3.70E+00
-3.72E+00
-3.74E+00
-3.75E+00
-3.77E+00
-3.79E+00
-3.80E+00
-3.82E+00
-3.84E+00
-3.86E+00
-3.87E+00
-3.89E+00
-3.91E+00
-3.93E+00
-3.94E+00
-3.96E+00
-3.98E+00
-4.00E+00
-4.01E+00
-4.03E+00
-4.05E+00
-4.07E+00
-4.08E+00
-4.10E+00
-4.12E+00
];

ADAMS_T_FRIC_1=-10;

ADAMS_M_FRIC_1=[8.50E+00
8.38E+00
8.16E+00
7.95E+00
7.74E+00
7.53E+00
7.33E+00
7.14E+00
6.95E+00
6.76E+00
6.57E+00
6.39E+00
6.21E+00
6.03E+00
5.86E+00
5.68E+00
5.51E+00
5.34E+00
5.18E+00
5.01E+00
4.85E+00
4.68E+00
4.52E+00
4.36E+00
4.21E+00
4.05E+00
3.89E+00
3.74E+00
3.58E+00
3.43E+00
3.27E+00
3.12E+00
2.97E+00
2.82E+00
2.67E+00
2.52E+00
2.37E+00
2.22E+00
2.07E+00
1.92E+00
1.77E+00
1.62E+00
1.47E+00
1.32E+00
1.17E+00
1.02E+00
8.70E-01
7.20E-01
5.69E-01
4.18E-01
2.67E-01
1.16E-01
-3.60E-02
-1.88E-01
-3.41E-01
-4.94E-01
-6.47E-01
-8.01E-01
-9.55E-01
-1.11E+00
-1.27E+00
-1.49E+00
-1.65E+00
-1.81E+00
-1.97E+00
-2.14E+00
-2.30E+00
-2.46E+00
-2.63E+00
-2.79E+00
-2.95E+00
-3.12E+00
-3.28E+00
-3.44E+00
-3.61E+00
-3.77E+00
-3.93E+00
-4.09E+00
-4.25E+00
-4.41E+00
-4.57E+00
-4.73E+00
-4.89E+00
-5.05E+00
-5.20E+00
-5.36E+00
-5.51E+00
-5.66E+00
-5.81E+00
-5.96E+00
-6.11E+00
-6.25E+00
-6.39E+00
-6.53E+00
-6.67E+00
-6.80E+00
-6.94E+00
-7.07E+00
-7.20E+00
-7.32E+00
-7.45E+00
-7.57E+00
-7.69E+00
-7.81E+00
-7.93E+00
-8.04E+00
-8.15E+00
-8.25E+00
-8.35E+00
-8.45E+00
-8.53E+00
-8.62E+00
-8.70E+00
-8.77E+00
-8.84E+00
-8.91E+00
-8.97E+00
-9.03E+00
-9.08E+00
-9.13E+00
-9.18E+00
-9.22E+00
-9.26E+00
-9.29E+00
-9.32E+00
-9.34E+00
-9.36E+00
-9.38E+00
-9.39E+00
-9.40E+00
-9.40E+00
-9.40E+00
-9.40E+00
-9.40E+00
-9.39E+00
-9.38E+00
-9.37E+00
-9.35E+00
-9.34E+00
-9.32E+00
-9.29E+00
-9.27E+00
-9.24E+00
-9.22E+00
-9.19E+00
-9.16E+00
-9.13E+00
-9.09E+00
-9.06E+00
-9.03E+00
-8.99E+00
-8.95E+00
-8.92E+00
-8.88E+00
-8.84E+00
-8.81E+00
-8.77E+00
-8.73E+00
-8.69E+00
-8.66E+00
-8.62E+00
-8.58E+00
-8.55E+00
-8.51E+00
-8.48E+00
-8.44E+00
-8.41E+00
-8.38E+00
-8.35E+00
-8.32E+00
-8.29E+00
-8.26E+00
-8.23E+00
-8.21E+00
-8.18E+00
-8.16E+00
-8.14E+00
-8.12E+00
-8.10E+00
-8.09E+00
-8.07E+00
-8.06E+00
-8.05E+00
-8.04E+00
-8.03E+00
-8.03E+00
-8.03E+00
-8.03E+00
-8.04E+00
-8.04E+00
-8.05E+00
-8.07E+00
-8.08E+00
-8.10E+00
-8.12E+00
-8.15E+00
-8.18E+00
-8.22E+00
-8.26E+00
-8.30E+00
-8.35E+00
-8.41E+00
-8.47E+00
-8.53E+00
-8.61E+00
-8.69E+00
-8.78E+00
-8.87E+00
-8.98E+00
-9.10E+00
-9.23E+00
-9.37E+00
-9.52E+00
-9.69E+00
-9.88E+00
-1.01E+01
-1.03E+01
-1.06E+01
-1.09E+01
-1.12E+01
-1.15E+01
-1.19E+01
-1.24E+01
-1.29E+01
-1.35E+01
-1.43E+01
-1.51E+01
-1.61E+01
-1.74E+01
-1.89E+01
-2.08E+01
-2.33E+01
-2.68E+01
-3.20E+01
-4.07E+01
-5.83E+01
-1.17E+02
];

ADAMS_a_FRIC_2=[0.00E+00
-1.75E-02
-3.49E-02
-5.24E-02
-6.98E-02
-8.73E-02
-1.05E-01
-1.22E-01
-1.40E-01
-1.57E-01
-1.75E-01
-1.92E-01
-2.09E-01
-2.27E-01
-2.44E-01
-2.62E-01
-2.79E-01
-2.97E-01
-3.14E-01
-3.32E-01
-3.49E-01
-3.67E-01
-3.84E-01
-4.01E-01
-4.19E-01
-4.36E-01
-4.54E-01
-4.71E-01
-4.89E-01
-5.06E-01
-5.24E-01
-5.41E-01
-5.59E-01
-5.76E-01
-5.93E-01
-6.11E-01
-6.28E-01
-6.46E-01
-6.63E-01
-6.81E-01
-6.98E-01
-7.16E-01
-7.33E-01
-7.50E-01
-7.68E-01
-7.85E-01
-8.03E-01
-8.20E-01
-8.38E-01
-8.55E-01
-8.73E-01
-8.90E-01
-9.08E-01
-9.25E-01
-9.42E-01
-9.60E-01
-9.77E-01
-9.95E-01
-1.01E+00
-1.03E+00
-1.05E+00
-1.06E+00
-1.08E+00
-1.10E+00
-1.12E+00
-1.13E+00
-1.15E+00
-1.17E+00
-1.19E+00
-1.20E+00
-1.22E+00
-1.24E+00
-1.26E+00
-1.27E+00
-1.29E+00
-1.31E+00
-1.33E+00
-1.34E+00
-1.36E+00
-1.38E+00
-1.40E+00
-1.41E+00
-1.43E+00
-1.45E+00
-1.47E+00
-1.48E+00
-1.50E+00
-1.52E+00
-1.54E+00
-1.55E+00
-1.57E+00
-1.59E+00
-1.61E+00
-1.62E+00
-1.64E+00
-1.66E+00
-1.68E+00
-1.69E+00
-1.71E+00
-1.73E+00
-1.75E+00
-1.76E+00
-1.78E+00
-1.80E+00
-1.82E+00
-1.83E+00
-1.85E+00
-1.87E+00
-1.88E+00
-1.90E+00
-1.92E+00
-1.94E+00
-1.95E+00
-1.97E+00
-1.99E+00
-2.01E+00
-2.02E+00
-2.04E+00
-2.06E+00
-2.08E+00
-2.09E+00
-2.11E+00
-2.13E+00
-2.15E+00
-2.16E+00
-2.18E+00
-2.20E+00
-2.22E+00
-2.23E+00
-2.25E+00
-2.27E+00
-2.29E+00
-2.30E+00
-2.32E+00
-2.34E+00
-2.36E+00
-2.37E+00
-2.39E+00
-2.41E+00
-2.43E+00
-2.44E+00
-2.46E+00
-2.48E+00
-2.50E+00
-2.51E+00
-2.53E+00
-2.55E+00
-2.57E+00
-2.58E+00
-2.60E+00
-2.62E+00
-2.64E+00
-2.65E+00
-2.67E+00
-2.69E+00
-2.71E+00
-2.72E+00
-2.74E+00
-2.76E+00
-2.78E+00
-2.79E+00
-2.81E+00
-2.83E+00
-2.84E+00
-2.86E+00
-2.88E+00
-2.90E+00
-2.91E+00
-2.93E+00
-2.95E+00
-2.97E+00
-2.98E+00
-3.00E+00
-3.02E+00
-3.04E+00
-3.05E+00
-3.07E+00
-3.09E+00
-3.11E+00
-3.12E+00
-3.14E+00
-3.16E+00
-3.18E+00
-3.19E+00
-3.21E+00
-3.23E+00
-3.25E+00
-3.26E+00
-3.28E+00
-3.30E+00
-3.32E+00
-3.33E+00
-3.35E+00
-3.37E+00
-3.39E+00
-3.40E+00
-3.42E+00
-3.44E+00
-3.46E+00
-3.47E+00
-3.49E+00
-3.51E+00
-3.53E+00
-3.54E+00
-3.56E+00
-3.58E+00
-3.60E+00
-3.61E+00
-3.63E+00
-3.65E+00
-3.67E+00
-3.68E+00
-3.70E+00
-3.72E+00
-3.74E+00
-3.75E+00
-3.77E+00
-3.79E+00
-3.80E+00
-3.82E+00
-3.84E+00
-3.86E+00
-3.87E+00
-3.89E+00
-3.91E+00
-3.93E+00
-3.94E+00
-3.96E+00
-3.98E+00
-4.00E+00
-4.01E+00
-4.03E+00
-4.05E+00
-4.07E+00
-4.08E+00
];

ADAMS_T_FRIC_2=-10;

ADAMS_M_FRIC_2=[7.60E+00
7.48E+00
7.31E+00
7.13E+00
6.96E+00
6.79E+00
6.62E+00
6.46E+00
6.30E+00
6.14E+00
5.98E+00
5.82E+00
5.67E+00
5.52E+00
5.37E+00
5.21E+00
5.07E+00
4.92E+00
4.77E+00
4.62E+00
4.48E+00
4.34E+00
4.19E+00
4.05E+00
3.91E+00
3.76E+00
3.62E+00
3.48E+00
3.34E+00
3.20E+00
3.06E+00
2.92E+00
2.78E+00
2.64E+00
2.50E+00
2.37E+00
2.23E+00
2.09E+00
1.95E+00
1.81E+00
1.67E+00
1.53E+00
1.39E+00
1.25E+00
1.11E+00
9.67E-01
8.25E-01
6.83E-01
5.41E-01
3.98E-01
2.54E-01
1.10E-01
-3.43E-02
-1.79E-01
-3.25E-01
-4.71E-01
-6.18E-01
-7.66E-01
-9.14E-01
-1.06E+00
-1.21E+00
-1.56E+00
-1.73E+00
-1.90E+00
-2.06E+00
-2.23E+00
-2.40E+00
-2.57E+00
-2.74E+00
-2.91E+00
-3.08E+00
-3.25E+00
-3.42E+00
-3.59E+00
-3.75E+00
-3.92E+00
-4.09E+00
-4.26E+00
-4.42E+00
-4.59E+00
-4.75E+00
-4.92E+00
-5.08E+00
-5.24E+00
-5.41E+00
-5.57E+00
-5.73E+00
-5.88E+00
-6.04E+00
-6.19E+00
-6.35E+00
-6.50E+00
-6.65E+00
-6.79E+00
-6.94E+00
-7.08E+00
-7.22E+00
-7.36E+00
-7.50E+00
-7.63E+00
-7.76E+00
-7.89E+00
-8.02E+00
-8.14E+00
-8.27E+00
-8.38E+00
-8.49E+00
-8.59E+00
-8.68E+00
-8.77E+00
-8.85E+00
-8.93E+00
-9.01E+00
-9.08E+00
-9.14E+00
-9.21E+00
-9.27E+00
-9.32E+00
-9.37E+00
-9.42E+00
-9.46E+00
-9.50E+00
-9.54E+00
-9.57E+00
-9.60E+00
-9.62E+00
-9.64E+00
-9.65E+00
-9.66E+00
-9.67E+00
-9.67E+00
-9.68E+00
-9.68E+00
-9.67E+00
-9.66E+00
-9.66E+00
-9.64E+00
-9.63E+00
-9.61E+00
-9.59E+00
-9.57E+00
-9.55E+00
-9.53E+00
-9.50E+00
-9.48E+00
-9.45E+00
-9.42E+00
-9.39E+00
-9.36E+00
-9.33E+00
-9.30E+00
-9.26E+00
-9.23E+00
-9.20E+00
-9.16E+00
-9.13E+00
-9.10E+00
-9.07E+00
-9.03E+00
-9.00E+00
-8.97E+00
-8.94E+00
-8.91E+00
-8.88E+00
-8.85E+00
-8.82E+00
-8.80E+00
-8.77E+00
-8.75E+00
-8.72E+00
-8.70E+00
-8.68E+00
-8.66E+00
-8.64E+00
-8.63E+00
-8.61E+00
-8.60E+00
-8.59E+00
-8.58E+00
-8.58E+00
-8.57E+00
-8.57E+00
-8.57E+00
-8.58E+00
-8.58E+00
-8.59E+00
-8.61E+00
-8.62E+00
-8.64E+00
-8.66E+00
-8.69E+00
-8.72E+00
-8.75E+00
-8.79E+00
-8.84E+00
-8.88E+00
-8.94E+00
-9.00E+00
-9.06E+00
-9.13E+00
-9.21E+00
-9.30E+00
-9.39E+00
-9.50E+00
-9.61E+00
-9.73E+00
-9.86E+00
-1.00E+01
-1.02E+01
-1.03E+01
-1.05E+01
-1.07E+01
-1.10E+01
-1.12E+01
-1.15E+01
-1.18E+01
-1.22E+01
-1.26E+01
-1.30E+01
-1.35E+01
-1.41E+01
-1.47E+01
-1.55E+01
-1.64E+01
-1.75E+01
-1.88E+01
-2.04E+01
-2.24E+01
-2.50E+01
-2.84E+01
-3.31E+01
-4.01E+01
-5.20E+01
-7.61E+01
-1.54E+02
];

ADAMS_a_FRIC_3=[0.00E+00
-1.75E-02
-3.49E-02
-5.24E-02
-6.98E-02
-8.73E-02
-1.05E-01
-1.22E-01
-1.40E-01
-1.57E-01
-1.75E-01
-1.92E-01
-2.09E-01
-2.27E-01
-2.44E-01
-2.62E-01
-2.79E-01
-2.97E-01
-3.14E-01
-3.32E-01
-3.49E-01
-3.67E-01
-3.84E-01
-4.01E-01
-4.19E-01
-4.36E-01
-4.54E-01
-4.71E-01
-4.89E-01
-5.06E-01
-5.24E-01
-5.41E-01
-5.59E-01
-5.76E-01
-5.93E-01
-6.11E-01
-6.28E-01
-6.46E-01
-6.63E-01
-6.81E-01
-6.98E-01
-7.16E-01
-7.33E-01
-7.50E-01
-7.68E-01
-7.85E-01
-8.03E-01
-8.20E-01
-8.38E-01
-8.55E-01
-8.73E-01
-8.90E-01
-9.08E-01
-9.25E-01
-9.42E-01
-9.60E-01
-9.77E-01
-9.95E-01
-1.01E+00
-1.03E+00
-1.05E+00
-1.06E+00
-1.08E+00
-1.10E+00
-1.12E+00
-1.13E+00
-1.15E+00
-1.17E+00
-1.19E+00
-1.20E+00
-1.22E+00
-1.24E+00
-1.26E+00
-1.27E+00
-1.29E+00
-1.31E+00
-1.33E+00
-1.34E+00
-1.36E+00
-1.38E+00
-1.40E+00
-1.41E+00
-1.43E+00
-1.45E+00
-1.47E+00
-1.48E+00
-1.50E+00
-1.52E+00
-1.54E+00
-1.55E+00
-1.57E+00
-1.59E+00
-1.61E+00
-1.62E+00
-1.64E+00
-1.66E+00
-1.68E+00
-1.69E+00
-1.71E+00
-1.73E+00
-1.75E+00
-1.76E+00
-1.78E+00
-1.80E+00
-1.82E+00
-1.83E+00
-1.85E+00
-1.87E+00
-1.88E+00
-1.90E+00
-1.92E+00
-1.94E+00
-1.95E+00
-1.97E+00
-1.99E+00
-2.01E+00
-2.02E+00
-2.04E+00
-2.06E+00
-2.08E+00
-2.09E+00
-2.11E+00
-2.13E+00
-2.15E+00
-2.16E+00
-2.18E+00
-2.20E+00
-2.22E+00
-2.23E+00
-2.25E+00
-2.27E+00
-2.29E+00
-2.30E+00
-2.32E+00
-2.34E+00
-2.36E+00
-2.37E+00
-2.39E+00
-2.41E+00
-2.43E+00
-2.44E+00
-2.46E+00
-2.48E+00
-2.50E+00
-2.51E+00
-2.53E+00
-2.55E+00
-2.57E+00
-2.58E+00
-2.60E+00
-2.62E+00
-2.64E+00
-2.65E+00
-2.67E+00
-2.69E+00
-2.71E+00
-2.72E+00
-2.74E+00
-2.76E+00
-2.78E+00
-2.79E+00
-2.81E+00
-2.83E+00
-2.84E+00
-2.86E+00
-2.88E+00
-2.90E+00
-2.91E+00
-2.93E+00
-2.95E+00
-2.97E+00
-2.98E+00
-3.00E+00
-3.02E+00
-3.04E+00
-3.05E+00
-3.07E+00
-3.09E+00
-3.11E+00
-3.12E+00
-3.14E+00
-3.16E+00
-3.18E+00
-3.19E+00
-3.21E+00
-3.23E+00
-3.25E+00
-3.26E+00
-3.28E+00
-3.30E+00
-3.32E+00
-3.33E+00
-3.35E+00
-3.37E+00
-3.39E+00
-3.40E+00
-3.42E+00
-3.44E+00
-3.46E+00
-3.47E+00
-3.49E+00
-3.51E+00
-3.53E+00
-3.54E+00
-3.56E+00
-3.58E+00
-3.60E+00
-3.61E+00
-3.63E+00
-3.65E+00
-3.67E+00
-3.68E+00
-3.70E+00
-3.72E+00
-3.74E+00
-3.75E+00
-3.77E+00
-3.79E+00
-3.80E+00
-3.82E+00
-3.84E+00
-3.86E+00
-3.87E+00
-3.89E+00
-3.91E+00
-3.93E+00
-3.94E+00
-3.96E+00
-3.98E+00
-4.00E+00
-4.01E+00
-4.03E+00
-4.05E+00
-4.07E+00
-4.08E+00
-4.10E+00
-4.12E+00
];

ADAMS_T_FRIC_3=-10;

ADAMS_M_FRIC_3=[7.15E+00
7.01E+00
6.80E+00
6.59E+00
6.39E+00
6.19E+00
6.00E+00
5.81E+00
5.62E+00
5.44E+00
5.26E+00
5.08E+00
4.91E+00
4.73E+00
4.56E+00
4.39E+00
4.22E+00
4.06E+00
3.89E+00
3.73E+00
3.57E+00
3.41E+00
3.25E+00
3.09E+00
2.93E+00
2.78E+00
2.62E+00
2.47E+00
2.31E+00
2.16E+00
2.01E+00
1.86E+00
1.70E+00
1.55E+00
1.40E+00
1.25E+00
1.10E+00
9.50E-01
8.00E-01
6.50E-01
5.00E-01
3.49E-01
1.99E-01
4.88E-02
-1.02E-01
-2.53E-01
-4.04E-01
-5.55E-01
-7.07E-01
-8.59E-01
-1.01E+00
-1.16E+00
-1.32E+00
-1.47E+00
-1.63E+00
-1.78E+00
-1.93E+00
-2.09E+00
-2.25E+00
-2.40E+00
-2.56E+00
-2.84E+00
-3.01E+00
-3.17E+00
-3.33E+00
-3.50E+00
-3.66E+00
-3.83E+00
-3.99E+00
-4.15E+00
-4.32E+00
-4.48E+00
-4.65E+00
-4.81E+00
-4.98E+00
-5.14E+00
-5.30E+00
-5.46E+00
-5.63E+00
-5.79E+00
-5.95E+00
-6.11E+00
-6.26E+00
-6.42E+00
-6.58E+00
-6.73E+00
-6.88E+00
-7.04E+00
-7.18E+00
-7.33E+00
-7.48E+00
-7.62E+00
-7.76E+00
-7.90E+00
-8.04E+00
-8.18E+00
-8.31E+00
-8.44E+00
-8.57E+00
-8.70E+00
-8.83E+00
-8.96E+00
-9.09E+00
-9.22E+00
-9.34E+00
-9.47E+00
-9.58E+00
-9.69E+00
-9.79E+00
-9.88E+00
-9.97E+00
-1.00E+01
-1.01E+01
-1.02E+01
-1.02E+01
-1.03E+01
-1.03E+01
-1.04E+01
-1.04E+01
-1.05E+01
-1.05E+01
-1.05E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.05E+01
-1.05E+01
-1.05E+01
-1.05E+01
-1.04E+01
-1.04E+01
-1.04E+01
-1.03E+01
-1.03E+01
-1.03E+01
-1.02E+01
-1.02E+01
-1.01E+01
-1.01E+01
-1.01E+01
-1.00E+01
-9.98E+00
-9.93E+00
-9.89E+00
-9.85E+00
-9.80E+00
-9.76E+00
-9.72E+00
-9.67E+00
-9.63E+00
-9.59E+00
-9.55E+00
-9.51E+00
-9.47E+00
-9.43E+00
-9.39E+00
-9.36E+00
-9.32E+00
-9.29E+00
-9.25E+00
-9.22E+00
-9.19E+00
-9.16E+00
-9.14E+00
-9.11E+00
-9.09E+00
-9.06E+00
-9.04E+00
-9.03E+00
-9.01E+00
-9.00E+00
-8.99E+00
-8.98E+00
-8.97E+00
-8.97E+00
-8.97E+00
-8.97E+00
-8.97E+00
-8.98E+00
-8.99E+00
-9.01E+00
-9.03E+00
-9.05E+00
-9.07E+00
-9.10E+00
-9.14E+00
-9.18E+00
-9.22E+00
-9.27E+00
-9.32E+00
-9.38E+00
-9.45E+00
-9.52E+00
-9.60E+00
-9.69E+00
-9.79E+00
-9.89E+00
-1.00E+01
-1.01E+01
-1.03E+01
-1.04E+01
-1.06E+01
-1.08E+01
-1.10E+01
-1.12E+01
-1.14E+01
-1.17E+01
-1.20E+01
-1.23E+01
-1.27E+01
-1.31E+01
-1.36E+01
-1.41E+01
-1.47E+01
-1.54E+01
-1.63E+01
-1.72E+01
-1.84E+01
-1.98E+01
-2.16E+01
-2.38E+01
-2.67E+01
-3.07E+01
-3.66E+01
-4.65E+01
-6.67E+01
-1.34E+02
];

ADAMS_a_FRIC_4=[0.00E+00
-1.75E-02
-3.49E-02
-5.24E-02
-6.98E-02
-8.73E-02
-1.05E-01
-1.22E-01
-1.40E-01
-1.57E-01
-1.75E-01
-1.92E-01
-2.09E-01
-2.27E-01
-2.44E-01
-2.62E-01
-2.79E-01
-2.97E-01
-3.14E-01
-3.32E-01
-3.49E-01
-3.67E-01
-3.84E-01
-4.01E-01
-4.19E-01
-4.36E-01
-4.54E-01
-4.71E-01
-4.89E-01
-5.06E-01
-5.24E-01
-5.41E-01
-5.59E-01
-5.76E-01
-5.93E-01
-6.11E-01
-6.28E-01
-6.46E-01
-6.63E-01
-6.81E-01
-6.98E-01
-7.16E-01
-7.33E-01
-7.50E-01
-7.68E-01
-7.85E-01
-8.03E-01
-8.20E-01
-8.38E-01
-8.55E-01
-8.73E-01
-8.90E-01
-9.08E-01
-9.25E-01
-9.42E-01
-9.60E-01
-9.77E-01
-9.95E-01
-1.01E+00
-1.03E+00
-1.05E+00
-1.06E+00
-1.08E+00
-1.10E+00
-1.12E+00
-1.13E+00
-1.15E+00
-1.17E+00
-1.19E+00
-1.20E+00
-1.22E+00
-1.24E+00
-1.26E+00
-1.27E+00
-1.29E+00
-1.31E+00
-1.33E+00
-1.34E+00
-1.36E+00
-1.38E+00
-1.40E+00
-1.41E+00
-1.43E+00
-1.45E+00
-1.47E+00
-1.48E+00
-1.50E+00
-1.52E+00
-1.54E+00
-1.55E+00
-1.57E+00
-1.59E+00
-1.61E+00
-1.62E+00
-1.64E+00
-1.66E+00
-1.68E+00
-1.69E+00
-1.71E+00
-1.73E+00
-1.75E+00
-1.76E+00
-1.78E+00
-1.80E+00
-1.82E+00
-1.83E+00
-1.85E+00
-1.87E+00
-1.88E+00
-1.90E+00
-1.92E+00
-1.94E+00
-1.95E+00
-1.97E+00
-1.99E+00
-2.01E+00
-2.02E+00
-2.04E+00
-2.06E+00
-2.08E+00
-2.09E+00
-2.11E+00
-2.13E+00
-2.15E+00
-2.16E+00
-2.18E+00
-2.20E+00
-2.22E+00
-2.23E+00
-2.25E+00
-2.27E+00
-2.29E+00
-2.30E+00
-2.32E+00
-2.34E+00
-2.36E+00
-2.37E+00
-2.39E+00
-2.41E+00
-2.43E+00
-2.44E+00
-2.46E+00
-2.48E+00
-2.50E+00
-2.51E+00
-2.53E+00
-2.55E+00
-2.57E+00
-2.58E+00
-2.60E+00
-2.62E+00
-2.64E+00
-2.65E+00
-2.67E+00
-2.69E+00
-2.71E+00
-2.72E+00
-2.74E+00
-2.76E+00
-2.78E+00
-2.79E+00
-2.81E+00
-2.83E+00
-2.84E+00
-2.86E+00
-2.88E+00
-2.90E+00
-2.91E+00
-2.93E+00
-2.95E+00
-2.97E+00
-2.98E+00
-3.00E+00
-3.02E+00
-3.04E+00
-3.05E+00
-3.07E+00
-3.09E+00
-3.11E+00
-3.12E+00
-3.14E+00
-3.16E+00
-3.18E+00
-3.19E+00
-3.21E+00
-3.23E+00
-3.25E+00
-3.26E+00
-3.28E+00
-3.30E+00
-3.32E+00
-3.33E+00
-3.35E+00
-3.37E+00
-3.39E+00
-3.40E+00
-3.42E+00
-3.44E+00
-3.46E+00
-3.47E+00
-3.49E+00
-3.51E+00
-3.53E+00
-3.54E+00
-3.56E+00
-3.58E+00
-3.60E+00
-3.61E+00
-3.63E+00
-3.65E+00
-3.67E+00
-3.68E+00
-3.70E+00
-3.72E+00
-3.74E+00
-3.75E+00
-3.77E+00
-3.79E+00
-3.80E+00
-3.82E+00
-3.84E+00
-3.86E+00
-3.87E+00
-3.89E+00
-3.91E+00
-3.93E+00
-3.94E+00
-3.96E+00
-3.98E+00
-4.00E+00
-4.01E+00
-4.03E+00
-4.05E+00
-4.07E+00
-4.08E+00
];

ADAMS_T_FRIC_4=-10;

ADAMS_M_FRIC_4=[6.40E+00
6.26E+00
6.08E+00
5.91E+00
5.75E+00
5.58E+00
5.42E+00
5.26E+00
5.10E+00
4.94E+00
4.79E+00
4.63E+00
4.48E+00
4.33E+00
4.18E+00
4.03E+00
3.88E+00
3.73E+00
3.59E+00
3.44E+00
3.30E+00
3.15E+00
3.01E+00
2.87E+00
2.72E+00
2.58E+00
2.44E+00
2.30E+00
2.16E+00
2.02E+00
1.88E+00
1.74E+00
1.60E+00
1.46E+00
1.32E+00
1.18E+00
1.04E+00
8.95E-01
7.54E-01
6.13E-01
4.72E-01
3.30E-01
1.88E-01
4.62E-02
-9.64E-02
-2.39E-01
-3.83E-01
-5.27E-01
-6.71E-01
-8.16E-01
-9.62E-01
-1.11E+00
-1.26E+00
-1.40E+00
-1.55E+00
-1.70E+00
-1.85E+00
-2.00E+00
-2.15E+00
-2.30E+00
-2.45E+00
-2.98E+00
-3.15E+00
-3.32E+00
-3.49E+00
-3.66E+00
-3.83E+00
-4.00E+00
-4.17E+00
-4.33E+00
-4.50E+00
-4.67E+00
-4.84E+00
-5.01E+00
-5.18E+00
-5.35E+00
-5.52E+00
-5.68E+00
-5.85E+00
-6.02E+00
-6.18E+00
-6.34E+00
-6.51E+00
-6.67E+00
-6.83E+00
-6.99E+00
-7.15E+00
-7.31E+00
-7.47E+00
-7.62E+00
-7.77E+00
-7.92E+00
-8.07E+00
-8.22E+00
-8.36E+00
-8.51E+00
-8.65E+00
-8.79E+00
-8.93E+00
-9.07E+00
-9.20E+00
-9.34E+00
-9.48E+00
-9.61E+00
-9.74E+00
-9.87E+00
-9.98E+00
-1.01E+01
-1.02E+01
-1.03E+01
-1.03E+01
-1.04E+01
-1.05E+01
-1.05E+01
-1.06E+01
-1.06E+01
-1.07E+01
-1.07E+01
-1.08E+01
-1.08E+01
-1.08E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.09E+01
-1.08E+01
-1.08E+01
-1.08E+01
-1.08E+01
-1.07E+01
-1.07E+01
-1.07E+01
-1.06E+01
-1.06E+01
-1.06E+01
-1.05E+01
-1.05E+01
-1.04E+01
-1.04E+01
-1.04E+01
-1.03E+01
-1.03E+01
-1.02E+01
-1.02E+01
-1.02E+01
-1.01E+01
-1.01E+01
-1.01E+01
-1.00E+01
-9.98E+00
-9.94E+00
-9.91E+00
-9.88E+00
-9.84E+00
-9.81E+00
-9.78E+00
-9.76E+00
-9.73E+00
-9.70E+00
-9.68E+00
-9.66E+00
-9.64E+00
-9.62E+00
-9.61E+00
-9.59E+00
-9.58E+00
-9.57E+00
-9.56E+00
-9.56E+00
-9.56E+00
-9.56E+00
-9.56E+00
-9.57E+00
-9.58E+00
-9.60E+00
-9.61E+00
-9.63E+00
-9.66E+00
-9.69E+00
-9.72E+00
-9.76E+00
-9.80E+00
-9.85E+00
-9.90E+00
-9.96E+00
-1.00E+01
-1.01E+01
-1.02E+01
-1.03E+01
-1.04E+01
-1.05E+01
-1.06E+01
-1.07E+01
-1.08E+01
-1.10E+01
-1.11E+01
-1.13E+01
-1.15E+01
-1.17E+01
-1.19E+01
-1.21E+01
-1.24E+01
-1.27E+01
-1.30E+01
-1.34E+01
-1.38E+01
-1.43E+01
-1.48E+01
-1.53E+01
-1.60E+01
-1.68E+01
-1.77E+01
-1.87E+01
-1.99E+01
-2.15E+01
-2.33E+01
-2.56E+01
-2.86E+01
-3.24E+01
-3.78E+01
-4.59E+01
-5.94E+01
-8.70E+01
-1.76E+02
];

ADAMS_F_OVER_M_FRIC=[8.21E-03
8.35E-03
8.46E-03
8.57E-03
8.70E-03
8.84E-03
8.99E-03
9.15E-03
9.32E-03
9.50E-03
9.70E-03
9.90E-03
1.01E-02
1.04E-02
1.06E-02
1.09E-02
1.12E-02
1.15E-02
1.19E-02
1.22E-02
1.26E-02
1.30E-02
1.35E-02
1.40E-02
1.46E-02
1.52E-02
1.58E-02
1.65E-02
1.73E-02
1.82E-02
1.92E-02
2.04E-02
2.16E-02
2.31E-02
2.47E-02
2.66E-02
2.89E-02
3.16E-02
3.48E-02
3.88E-02
4.39E-02
5.06E-02
5.96E-02
7.27E-02
9.31E-02
1.30E-01
2.13E-01
6.04E-01
-7.22E-01
-2.26E-01
-1.34E-01
-9.49E-02
-7.36E-02
-6.01E-02
-5.08E-02
-4.39E-02
-3.87E-02
-3.46E-02
-3.13E-02
-2.86E-02
-2.63E-02
-2.43E-02
-2.27E-02
-2.12E-02
-1.99E-02
-1.88E-02
-1.78E-02
-1.69E-02
-1.61E-02
-1.53E-02
-1.47E-02
-1.41E-02
-1.35E-02
-1.30E-02
-1.25E-02
-1.21E-02
-1.17E-02
-1.13E-02
-1.10E-02
-1.06E-02
-1.03E-02
-1.01E-02
-9.79E-03
-9.54E-03
-9.31E-03
-9.09E-03
-8.89E-03
-8.69E-03
-8.51E-03
-8.34E-03
-8.18E-03
-8.03E-03
-7.89E-03
-7.76E-03
-7.64E-03
-7.52E-03
-7.42E-03
-7.32E-03
-7.23E-03
-7.15E-03
-7.08E-03
-7.02E-03
-6.97E-03
-6.92E-03
-6.88E-03
-6.85E-03
-6.81E-03
-6.77E-03
-6.73E-03
-6.68E-03
-6.62E-03
-6.56E-03
-6.50E-03
-6.44E-03
-6.39E-03
-6.33E-03
-6.28E-03
-6.23E-03
-6.18E-03
-6.14E-03
-6.10E-03
-6.06E-03
-6.03E-03
-6.00E-03
-5.97E-03
-5.94E-03
-5.92E-03
-5.89E-03
-5.87E-03
-5.86E-03
-5.84E-03
-5.83E-03
-5.82E-03
-5.81E-03
-5.80E-03
-5.80E-03
-5.79E-03
-5.79E-03
-5.79E-03
-5.79E-03
-5.79E-03
-5.80E-03
-5.80E-03
-5.81E-03
-5.82E-03
-5.83E-03
-5.84E-03
-5.85E-03
-5.86E-03
-5.88E-03
-5.89E-03
-5.91E-03
-5.93E-03
-5.95E-03
-5.97E-03
-5.99E-03
-6.01E-03
-6.03E-03
-6.06E-03
-6.08E-03
-6.11E-03
-6.13E-03
-6.16E-03
-6.19E-03
-6.22E-03
-6.25E-03
-6.28E-03
-6.31E-03
-6.34E-03
-6.37E-03
-6.40E-03
-6.44E-03
-6.47E-03
-6.51E-03
-6.54E-03
-6.58E-03
-6.62E-03
-6.65E-03
-6.69E-03
-6.73E-03
-6.77E-03
-6.81E-03
-6.85E-03
-6.89E-03
-6.93E-03
-6.97E-03
-7.01E-03
-7.05E-03
-7.09E-03
-7.13E-03
-7.18E-03
-7.22E-03
-7.26E-03
-7.30E-03
-7.35E-03
-7.39E-03
-7.43E-03
-7.48E-03
-7.52E-03
-7.56E-03
-7.61E-03
-7.65E-03
-7.69E-03
-7.74E-03
-7.78E-03
-7.82E-03
-7.86E-03
-7.91E-03
-7.95E-03
-7.99E-03
-8.03E-03
-8.07E-03
-8.11E-03
-8.15E-03
-8.19E-03
-8.23E-03
-8.26E-03
-8.30E-03
-8.34E-03
-8.38E-03
-8.41E-03
-8.45E-03
-8.49E-03
-8.53E-03
-8.57E-03
-8.61E-03
-8.65E-03
-8.67E-03
-8.69E-03
-8.69E-03
-8.68E-03
-8.66E-03
-8.63E-03
-8.59E-03
-8.54E-03
-8.46E-03
];

ADAMS_F_OVER_M_FRIC_1=[7.55E-03
7.66E-03
7.75E-03
7.84E-03
7.94E-03
8.05E-03
8.17E-03
8.30E-03
8.43E-03
8.57E-03
8.73E-03
8.89E-03
9.07E-03
9.25E-03
9.45E-03
9.66E-03
9.89E-03
1.01E-02
1.04E-02
1.07E-02
1.10E-02
1.13E-02
1.16E-02
1.20E-02
1.24E-02
1.28E-02
1.32E-02
1.37E-02
1.43E-02
1.49E-02
1.55E-02
1.62E-02
1.70E-02
1.79E-02
1.88E-02
1.99E-02
2.12E-02
2.26E-02
2.42E-02
2.60E-02
2.82E-02
3.08E-02
3.39E-02
3.77E-02
4.25E-02
4.87E-02
5.71E-02
6.90E-02
8.73E-02
1.19E-01
1.86E-01
4.30E-01
-1.38E+00
-2.65E-01
-1.46E-01
-1.01E-01
-7.73E-02
-6.25E-02
-5.25E-02
-4.52E-02
-3.97E-02
-3.54E-02
-3.20E-02
-2.92E-02
-2.68E-02
-2.48E-02
-2.31E-02
-2.16E-02
-2.03E-02
-1.91E-02
-1.81E-02
-1.72E-02
-1.63E-02
-1.56E-02
-1.49E-02
-1.43E-02
-1.37E-02
-1.32E-02
-1.27E-02
-1.23E-02
-1.19E-02
-1.15E-02
-1.12E-02
-1.08E-02
-1.05E-02
-1.03E-02
-9.99E-03
-9.75E-03
-9.52E-03
-9.30E-03
-9.10E-03
-8.92E-03
-8.74E-03
-8.58E-03
-8.42E-03
-8.28E-03
-8.15E-03
-8.03E-03
-7.92E-03
-7.82E-03
-7.74E-03
-7.66E-03
-7.60E-03
-7.54E-03
-7.50E-03
-7.46E-03
-7.41E-03
-7.37E-03
-7.31E-03
-7.25E-03
-7.18E-03
-7.10E-03
-7.03E-03
-6.96E-03
-6.89E-03
-6.82E-03
-6.76E-03
-6.70E-03
-6.64E-03
-6.59E-03
-6.54E-03
-6.50E-03
-6.45E-03
-6.41E-03
-6.38E-03
-6.34E-03
-6.31E-03
-6.29E-03
-6.26E-03
-6.24E-03
-6.22E-03
-6.20E-03
-6.19E-03
-6.17E-03
-6.16E-03
-6.15E-03
-6.15E-03
-6.14E-03
-6.14E-03
-6.14E-03
-6.14E-03
-6.14E-03
-6.15E-03
-6.15E-03
-6.16E-03
-6.17E-03
-6.18E-03
-6.19E-03
-6.20E-03
-6.21E-03
-6.23E-03
-6.25E-03
-6.27E-03
-6.28E-03
-6.30E-03
-6.33E-03
-6.35E-03
-6.37E-03
-6.40E-03
-6.42E-03
-6.45E-03
-6.48E-03
-6.51E-03
-6.54E-03
-6.57E-03
-6.60E-03
-6.63E-03
-6.66E-03
-6.70E-03
-6.73E-03
-6.77E-03
-6.80E-03
-6.84E-03
-6.88E-03
-6.91E-03
-6.95E-03
-6.99E-03
-7.03E-03
-7.07E-03
-7.11E-03
-7.16E-03
-7.20E-03
-7.24E-03
-7.29E-03
-7.33E-03
-7.37E-03
-7.42E-03
-7.46E-03
-7.51E-03
-7.56E-03
-7.60E-03
-7.65E-03
-7.69E-03
-7.74E-03
-7.79E-03
-7.84E-03
-7.88E-03
-7.93E-03
-7.98E-03
-8.03E-03
-8.08E-03
-8.12E-03
-8.17E-03
-8.22E-03
-8.27E-03
-8.32E-03
-8.36E-03
-8.41E-03
-8.46E-03
-8.50E-03
-8.55E-03
-8.60E-03
-8.64E-03
-8.69E-03
-8.73E-03
-8.77E-03
-8.82E-03
-8.86E-03
-8.90E-03
-8.94E-03
-8.99E-03
-9.03E-03
-9.07E-03
-9.12E-03
-9.17E-03
-9.22E-03
-9.26E-03
-9.29E-03
-9.30E-03
-9.31E-03
-9.30E-03
-9.28E-03
-9.25E-03
-9.21E-03
-9.15E-03
-9.07E-03
-8.96E-03
];

ADAMS_F_OVER_M_FRIC_2=[7.55E-03
7.66E-03
7.75E-03
7.84E-03
7.94E-03
8.05E-03
8.17E-03
8.30E-03
8.43E-03
8.57E-03
8.73E-03
8.89E-03
9.07E-03
9.25E-03
9.45E-03
9.66E-03
9.89E-03
1.01E-02
1.04E-02
1.07E-02
1.10E-02
1.13E-02
1.16E-02
1.20E-02
1.24E-02
1.28E-02
1.32E-02
1.37E-02
1.43E-02
1.49E-02
1.55E-02
1.62E-02
1.70E-02
1.79E-02
1.88E-02
1.99E-02
2.12E-02
2.26E-02
2.42E-02
2.60E-02
2.82E-02
3.08E-02
3.39E-02
3.77E-02
4.25E-02
4.87E-02
5.71E-02
6.90E-02
8.73E-02
1.19E-01
1.86E-01
4.30E-01
-1.38E+00
-2.65E-01
-1.46E-01
-1.01E-01
-7.73E-02
-6.25E-02
-5.25E-02
-4.52E-02
-3.97E-02
-3.54E-02
-3.20E-02
-2.92E-02
-2.68E-02
-2.48E-02
-2.31E-02
-2.16E-02
-2.03E-02
-1.91E-02
-1.81E-02
-1.72E-02
-1.63E-02
-1.56E-02
-1.49E-02
-1.43E-02
-1.37E-02
-1.32E-02
-1.27E-02
-1.23E-02
-1.19E-02
-1.15E-02
-1.12E-02
-1.08E-02
-1.05E-02
-1.03E-02
-9.99E-03
-9.75E-03
-9.52E-03
-9.30E-03
-9.10E-03
-8.92E-03
-8.74E-03
-8.58E-03
-8.42E-03
-8.28E-03
-8.15E-03
-8.03E-03
-7.92E-03
-7.82E-03
-7.74E-03
-7.66E-03
-7.60E-03
-7.54E-03
-7.50E-03
-7.46E-03
-7.41E-03
-7.37E-03
-7.31E-03
-7.25E-03
-7.18E-03
-7.10E-03
-7.03E-03
-6.96E-03
-6.89E-03
-6.82E-03
-6.76E-03
-6.70E-03
-6.64E-03
-6.59E-03
-6.54E-03
-6.50E-03
-6.45E-03
-6.41E-03
-6.38E-03
-6.34E-03
-6.31E-03
-6.29E-03
-6.26E-03
-6.24E-03
-6.22E-03
-6.20E-03
-6.19E-03
-6.17E-03
-6.16E-03
-6.15E-03
-6.15E-03
-6.14E-03
-6.14E-03
-6.14E-03
-6.14E-03
-6.14E-03
-6.15E-03
-6.15E-03
-6.16E-03
-6.17E-03
-6.18E-03
-6.19E-03
-6.20E-03
-6.21E-03
-6.23E-03
-6.25E-03
-6.27E-03
-6.28E-03
-6.30E-03
-6.33E-03
-6.35E-03
-6.37E-03
-6.40E-03
-6.42E-03
-6.45E-03
-6.48E-03
-6.51E-03
-6.54E-03
-6.57E-03
-6.60E-03
-6.63E-03
-6.66E-03
-6.70E-03
-6.73E-03
-6.77E-03
-6.80E-03
-6.84E-03
-6.88E-03
-6.91E-03
-6.95E-03
-6.99E-03
-7.03E-03
-7.07E-03
-7.11E-03
-7.16E-03
-7.20E-03
-7.24E-03
-7.29E-03
-7.33E-03
-7.37E-03
-7.42E-03
-7.46E-03
-7.51E-03
-7.56E-03
-7.60E-03
-7.65E-03
-7.69E-03
-7.74E-03
-7.79E-03
-7.84E-03
-7.88E-03
-7.93E-03
-7.98E-03
-8.03E-03
-8.08E-03
-8.12E-03
-8.17E-03
-8.22E-03
-8.27E-03
-8.32E-03
-8.36E-03
-8.41E-03
-8.46E-03
-8.50E-03
-8.55E-03
-8.60E-03
-8.64E-03
-8.69E-03
-8.73E-03
-8.77E-03
-8.82E-03
-8.86E-03
-8.90E-03
-8.94E-03
-8.99E-03
-9.03E-03
-9.07E-03
-9.12E-03
-9.17E-03
-9.22E-03
-9.26E-03
-9.29E-03
-9.30E-03
-9.31E-03
-9.30E-03
-9.28E-03
-9.25E-03
-9.21E-03
-9.15E-03
];

ADAMS_F_OVER_M_FRIC_3=[8.98E-03
9.16E-03
9.30E-03
9.46E-03
9.62E-03
9.80E-03
9.99E-03
1.02E-02
1.04E-02
1.06E-02
1.09E-02
1.12E-02
1.15E-02
1.18E-02
1.21E-02
1.25E-02
1.29E-02
1.33E-02
1.38E-02
1.43E-02
1.49E-02
1.55E-02
1.62E-02
1.69E-02
1.77E-02
1.86E-02
1.97E-02
2.08E-02
2.21E-02
2.36E-02
2.53E-02
2.73E-02
2.97E-02
3.25E-02
3.59E-02
4.01E-02
4.55E-02
5.26E-02
6.24E-02
7.67E-02
9.97E-02
1.42E-01
2.50E-01
1.02E+00
-4.88E-01
-1.97E-01
-1.23E-01
-8.95E-02
-7.03E-02
-5.79E-02
-4.92E-02
-4.27E-02
-3.78E-02
-3.39E-02
-3.07E-02
-2.81E-02
-2.58E-02
-2.39E-02
-2.23E-02
-2.09E-02
-1.96E-02
-1.85E-02
-1.75E-02
-1.67E-02
-1.59E-02
-1.51E-02
-1.45E-02
-1.39E-02
-1.33E-02
-1.28E-02
-1.24E-02
-1.19E-02
-1.15E-02
-1.11E-02
-1.08E-02
-1.05E-02
-1.02E-02
-9.89E-03
-9.63E-03
-9.38E-03
-9.14E-03
-8.92E-03
-8.72E-03
-8.52E-03
-8.34E-03
-8.16E-03
-8.00E-03
-7.85E-03
-7.70E-03
-7.56E-03
-7.43E-03
-7.31E-03
-7.20E-03
-7.09E-03
-6.99E-03
-6.89E-03
-6.81E-03
-6.73E-03
-6.65E-03
-6.59E-03
-6.53E-03
-6.47E-03
-6.43E-03
-6.39E-03
-6.36E-03
-6.33E-03
-6.30E-03
-6.27E-03
-6.23E-03
-6.19E-03
-6.14E-03
-6.09E-03
-6.05E-03
-6.00E-03
-5.95E-03
-5.90E-03
-5.86E-03
-5.82E-03
-5.78E-03
-5.75E-03
-5.71E-03
-5.68E-03
-5.66E-03
-5.63E-03
-5.61E-03
-5.58E-03
-5.57E-03
-5.55E-03
-5.53E-03
-5.52E-03
-5.51E-03
-5.50E-03
-5.49E-03
-5.48E-03
-5.48E-03
-5.48E-03
-5.47E-03
-5.47E-03
-5.48E-03
-5.48E-03
-5.48E-03
-5.49E-03
-5.49E-03
-5.50E-03
-5.51E-03
-5.52E-03
-5.53E-03
-5.55E-03
-5.56E-03
-5.57E-03
-5.59E-03
-5.61E-03
-5.62E-03
-5.64E-03
-5.66E-03
-5.68E-03
-5.70E-03
-5.73E-03
-5.75E-03
-5.77E-03
-5.80E-03
-5.82E-03
-5.85E-03
-5.88E-03
-5.90E-03
-5.93E-03
-5.96E-03
-5.99E-03
-6.02E-03
-6.05E-03
-6.08E-03
-6.11E-03
-6.14E-03
-6.18E-03
-6.21E-03
-6.24E-03
-6.28E-03
-6.31E-03
-6.35E-03
-6.38E-03
-6.42E-03
-6.46E-03
-6.49E-03
-6.53E-03
-6.57E-03
-6.60E-03
-6.64E-03
-6.68E-03
-6.72E-03
-6.76E-03
-6.80E-03
-6.83E-03
-6.87E-03
-6.91E-03
-6.95E-03
-6.99E-03
-7.03E-03
-7.07E-03
-7.11E-03
-7.15E-03
-7.19E-03
-7.23E-03
-7.27E-03
-7.31E-03
-7.34E-03
-7.38E-03
-7.42E-03
-7.46E-03
-7.50E-03
-7.53E-03
-7.57E-03
-7.61E-03
-7.64E-03
-7.68E-03
-7.71E-03
-7.75E-03
-7.78E-03
-7.81E-03
-7.84E-03
-7.88E-03
-7.91E-03
-7.94E-03
-7.98E-03
-8.01E-03
-8.05E-03
-8.08E-03
-8.11E-03
-8.14E-03
-8.15E-03
-8.15E-03
-8.14E-03
-8.12E-03
-8.09E-03
-8.05E-03
-8.00E-03
-7.93E-03
-7.83E-03
];

ADAMS_F_OVER_M_FRIC_4=[8.98E-03
9.16E-03
9.30E-03
9.46E-03
9.62E-03
9.80E-03
9.99E-03
1.02E-02
1.04E-02
1.06E-02
1.09E-02
1.12E-02
1.15E-02
1.18E-02
1.21E-02
1.25E-02
1.29E-02
1.33E-02
1.38E-02
1.43E-02
1.49E-02
1.55E-02
1.62E-02
1.69E-02
1.77E-02
1.86E-02
1.97E-02
2.08E-02
2.21E-02
2.36E-02
2.53E-02
2.73E-02
2.97E-02
3.25E-02
3.59E-02
4.01E-02
4.55E-02
5.26E-02
6.24E-02
7.67E-02
9.97E-02
1.42E-01
2.50E-01
1.02E+00
-4.88E-01
-1.97E-01
-1.23E-01
-8.95E-02
-7.03E-02
-5.79E-02
-4.92E-02
-4.27E-02
-3.78E-02
-3.39E-02
-3.07E-02
-2.81E-02
-2.58E-02
-2.39E-02
-2.23E-02
-2.09E-02
-1.96E-02
-1.85E-02
-1.75E-02
-1.67E-02
-1.59E-02
-1.51E-02
-1.45E-02
-1.39E-02
-1.33E-02
-1.28E-02
-1.24E-02
-1.19E-02
-1.15E-02
-1.11E-02
-1.08E-02
-1.05E-02
-1.02E-02
-9.89E-03
-9.63E-03
-9.38E-03
-9.14E-03
-8.92E-03
-8.72E-03
-8.52E-03
-8.34E-03
-8.16E-03
-8.00E-03
-7.85E-03
-7.70E-03
-7.56E-03
-7.43E-03
-7.31E-03
-7.20E-03
-7.09E-03
-6.99E-03
-6.89E-03
-6.81E-03
-6.73E-03
-6.65E-03
-6.59E-03
-6.53E-03
-6.47E-03
-6.43E-03
-6.39E-03
-6.36E-03
-6.33E-03
-6.30E-03
-6.27E-03
-6.23E-03
-6.19E-03
-6.14E-03
-6.09E-03
-6.05E-03
-6.00E-03
-5.95E-03
-5.90E-03
-5.86E-03
-5.82E-03
-5.78E-03
-5.75E-03
-5.71E-03
-5.68E-03
-5.66E-03
-5.63E-03
-5.61E-03
-5.58E-03
-5.57E-03
-5.55E-03
-5.53E-03
-5.52E-03
-5.51E-03
-5.50E-03
-5.49E-03
-5.48E-03
-5.48E-03
-5.48E-03
-5.47E-03
-5.47E-03
-5.48E-03
-5.48E-03
-5.48E-03
-5.49E-03
-5.49E-03
-5.50E-03
-5.51E-03
-5.52E-03
-5.53E-03
-5.55E-03
-5.56E-03
-5.57E-03
-5.59E-03
-5.61E-03
-5.62E-03
-5.64E-03
-5.66E-03
-5.68E-03
-5.70E-03
-5.73E-03
-5.75E-03
-5.77E-03
-5.80E-03
-5.82E-03
-5.85E-03
-5.88E-03
-5.90E-03
-5.93E-03
-5.96E-03
-5.99E-03
-6.02E-03
-6.05E-03
-6.08E-03
-6.11E-03
-6.14E-03
-6.18E-03
-6.21E-03
-6.24E-03
-6.28E-03
-6.31E-03
-6.35E-03
-6.38E-03
-6.42E-03
-6.46E-03
-6.49E-03
-6.53E-03
-6.57E-03
-6.60E-03
-6.64E-03
-6.68E-03
-6.72E-03
-6.76E-03
-6.80E-03
-6.83E-03
-6.87E-03
-6.91E-03
-6.95E-03
-6.99E-03
-7.03E-03
-7.07E-03
-7.11E-03
-7.15E-03
-7.19E-03
-7.23E-03
-7.27E-03
-7.31E-03
-7.34E-03
-7.38E-03
-7.42E-03
-7.46E-03
-7.50E-03
-7.53E-03
-7.57E-03
-7.61E-03
-7.64E-03
-7.68E-03
-7.71E-03
-7.75E-03
-7.78E-03
-7.81E-03
-7.84E-03
-7.88E-03
-7.91E-03
-7.94E-03
-7.98E-03
-8.01E-03
-8.05E-03
-8.08E-03
-8.11E-03
-8.14E-03
-8.15E-03
-8.15E-03
-8.14E-03
-8.12E-03
-8.09E-03
-8.05E-03
-8.00E-03
];
end