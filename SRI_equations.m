function [ LHS_box, RHS_box, Info_ ] = SRI_equations( Eqn_, EigFun_, T_, mu_, eta_, m_, k_, Re_, N_ )
%Dr Luke Robins 2019 luke.robins@cantab.net
%
%equations_ for GenEig for the viscous SRI, with no-slip Boundary
%Conditions.
%
%The equations represented here are those shown in the pdf
%"SRI_equations.pdf" in the viscous section.
%[L] = l = r_2-r_1 (lengthscale)
%[t] = 1 = 1/Omega_1 (timescale)
%
%Note that the term "co-location points" is used throughout this file. This
%refers to the Gauss-Lobatto quadrature points, albeit re-scaled to fit
%across the radial range r rather than across -1<x<+1. We use Gauss-Lobatto
%points ranging from 0 to T_ (i.e. x_j = cos(pi*j/T_).)
%
%Input Parameters
%  - Eqn_
%    This is an index corresponding to one of the equations described in
%    this file.
%      Eqn_=1 corresponds to the   radial Navier Stokes Equation (NS_r)
%      Eqn_=2 corresponds to the  angular Navier Stokes Equation (NS_theta)
%      Eqn_=3 corresponds to the vertical Navier Stokes Equation (NS_z)
%      Eqn_=4 corresponds to the mass equation.                  (Mass)
%      Eqn_=5 corresponds to the divergence-free equation.       (Div)
%
%   Boundary conditions are accessed using a negative equation index
%   Eqn_<0. The index EigFun_ is used as normal. Boundary conditions return
%   LHS_box and RHS_box as a horizontal vector, representing the evaluation
%   of terms on either of the boundaries (and at the corresponding
%   co-location point).
%       In this case there are 6 boundary conditions, from the requirement
%       that the three velocity perturbations go to zero on both
%       boundaries.
%
%  - EigFun_
%    This is an index corresponding to one of the Eigen-Functions that make
%    up each equation.
%      EigFun_=1 corresponds to the radial velocity perturbation u(r)
%      EigFun_=2 corresponds to the angular velocity perturbation v(r)
%      EigFun_=3 corresponds to the vertical velocity perturbation w(r)
%      EigFun_=4 corresponds to the density perturbation rho(r)
%      EigFun_=5 corresponds to the pressure perturbation P(r)
%
%  - See SRI_solver.m for a full breakdown the remaining input parameters.
%
%Output Parameters
%  - LHS_box
%    For a given equation Eqn_ and Eigen-Function EigFun_, LHS_box
%    is a matrix representing all terms involving Eigen-Function EigFun_
%    that are on the Left-Hand-Side of the equation Eqn_. This matrix is
%    not required to be square.
%       Specifically, LHS_box is a co-location matrix. It is set up to
%    multiply with a vector representing EigFun_, returning a vector
%    corresponding to values at each of the co-location points. These
%    values are the summed contributions of the LHS terms involving
%    Eigen-Function EigFun_ at each of the co-location points.
%           (For example, in the mass equation (Eqn_=4) for the 
%            Eigen-Function w (EigFun_=3) the values returned would be N_^2
%            times the value of the vertical velocity perturbation at each
%            co-location point.
%            On the other hand, in the divergence-free equation (Eqn_=5)
%            for the Eigen-Function u (EigFun_=1) the values returned would
%            be (du/dr + u/r) evaluated at each co-location point.
%
%  - RHS_box
%    The same as LHS_box, but for the Right-Hand-Side of each equation.
%    Note that the equations have been formulated such that any terms
%    multiplied by the eigenvalue (the complex growth rate sigma) are on
%    the right hand side, and all other terms are on the LHS.
%
%  - Info_
%    Generally Info_=0. However, it is sometimes necessary to pass on
%    information about the equations, rather than terms within them. For
%    example, the number of equations (which is also the number of
%    Eigen-Functions). This information can be retrieved by using specific
%    non-standard combinations of indices. 
%    There are four such cases: 
%    1) [~,~,Info_] = equations_(0,0,0);
%       Here Info_ will be the number of equations in the system; which is
%       necessarily also the number of Eigen-Functions in the system.
%
%    2) How many boundary conditions are in the system?
%       [~,~,Info_] = equations_(-1,0,0);
%       Here Info_ will be the number of boundary conditions within the
%       system.
%
%    3) Which equations should not be evaluated on the boundaries?
%           This may be necessary if an equation is made irrelevant at the
%           boundary by the boundary condition.
%           This is evaluated for each equation eqn_ and for each of the
%           two boundaries b_=+/-1. b_=-1 denotes the inner boundary, b_=+1
%           denotes the outer.
%       [~,~,Info_]=equations(eqn_,0,b_);
%       Here Info_=1 means you should not evaluate the equation
%       corresponding to eqn_ on the boundary. Otherwise Info_=0.
%
%    4) Which Eigen-Functions require variations on their number of
%       coefficients?
%       [~,~,Info_]=equations(0,EigFun_,0);
%       Here Info_ will be the modification to the number of Chebyshev
%       expansion coefficients that the Eigen-Function EigFun_ has. This is
%       refered to as C_mod in GenEig.m
%       C_mod can be positive or negative; the total is then:
%       1+Num_Cheb_terms+C_mod.
%       (The +1 handles the zeroth coefficient.)
%           These modifications exist to ensure that the larger matrix can
%           be kept square, despite any imposed boundary conditions adding
%           further lines.
%           These are typically unnecessary, hence C_mod is generally 0.
%
Z_=T_+1;
%Z_ is purely for book-keeping. T_ is the term at which we truncate the
%Chebyshev expansion. However, including the zeroth term of the expansion,
%this means that there are actually Z_=T_+1 terms in an expansion truncated
%at this point.

%The number of equations. This is also the number of Eigen-Functions.
N_Eqns_EigFuns=5;

%The number of boundary conditions.
N_BCs=6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialise:
LHS_box=0;
RHS_box=0;
Info_=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Co-location points are common for every equation:
cs_=zeros(Z_,1);
for a_=1:Z_
    j_=a_-1;
    %Co-location points:
    cs_(a_)=cos(j_*pi/T_);
end
%We rescale the co-location points to range from 0 (the inner boundary) to
%1 (the outer boundary): 
x_min=0;
x_max=1;
%Co-location points scaled across the gap-width:
xm_=-cs_*((x_max-x_min)/2)+((x_max+x_min)/2);
dcdx_ = -2/(x_max-x_min);
%The minus sign ensures that xm_(a_=1) corresponds to the INNER boundary,
%rather than outer.

%Radius:
r_=(eta_/(1-eta_))+xm_;

%Basic state angular velocity of a Taylor-Couette system:
A_=(mu_-eta_^2)/(1-(eta_^2));
B_=((eta_/(1-eta_))^2)*(1-mu_)/(1-(eta_^2));
Omega_=A_+B_./(r_.*r_);

%Both the radius and Omega are vectors, giving their respective quantities
%at each co-location point.

%The vorticity of this system is a constant:
zeta_=2*A_;

%Kinematic Viscosity:
%(Expressed using the length- and time-scales described above)
nu_=(1/Re_)*(eta_/(1-eta_));

%The size of each Eigen-Function's Chebyshev Expansion:
box_width=[Z_;Z_;Z_;Z_;Z_];
%This is returned as C_Mod in SRI_solver.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specific Info:
%For specific inputs, specific information can be retrieved from an
%equations_ function rather than its usual output:
if Eqn_==0
    if EigFun_==0
        %Number of equations/Eigen-Functions
        Info_=N_Eqns_EigFuns;
        return
    end
    if EigFun_>0
        %Return any modification to the Chebyshev Expansion
        Info_=box_width(EigFun_)-Z_;
        return
    end
    return
end
if and(Eqn_==-1,EigFun_==0)
    %Return the number of boundary conditions.
    Info_=N_BCs;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NS1:
%Radial Navier Stokes equation:
if Eqn_==1
    
    %u
    if EigFun_==1
        u_=CMatrix(Z_,box_width(EigFun_));
        Du_=u_*dcdx_*DMatrix(box_width(EigFun_));
        DDu_=Du_*dcdx_*DMatrix(box_width(EigFun_));
        
        LHS_box=nu_*DDu_+fT_((nu_./r_),Du_)-fT_(1i*m_*Omega_+nu_*(k_^2+((m_^2+1)./(r_.*r_))),u_);
        RHS_box=u_;
    end
    
    %v
    if EigFun_==2
        v_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=fT_(2*Omega_-(nu_*(2i*m_)./(r_.*r_)),v_);
    end
    
    %w
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %rho
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %P
    if EigFun_==5
        P_=CMatrix(Z_,box_width(EigFun_));
        DP_=P_*dcdx_*DMatrix(box_width(EigFun_));
        
        LHS_box=-DP_;
    end
    
    %Do not evaluate on the boundaries:
    if and(EigFun_==0,or(T_==1,T_==-1))
        Info_=1;
        return
    end
    %Trim boundaries:
    %We don't wish to evaluate this equation on the boundaries because of
    %the boundary conditions.
    if not(length(LHS_box)==1)
        LHS_box=LHS_box(2:Z_-1,:);
    end
    if not(length(RHS_box)==1)
        RHS_box=RHS_box(2:Z_-1,:);
    end
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NS2:
%Angular Navier Stokes Equation:
if Eqn_==2
    
    %u
    if EigFun_==1
        u_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=fT_(-zeta_+(2i*nu_*m_./(r_.*r_)),u_);
    end
    
    %v
    if EigFun_==2
        v_=CMatrix(Z_,box_width(EigFun_));
        Dv_=v_*dcdx_*DMatrix(box_width(EigFun_));
        DDv_=Dv_*dcdx_*DMatrix(box_width(EigFun_));
        
        LHS_box=nu_*DDv_+fT_((nu_./r_),Dv_)-fT_(1i*m_*Omega_+nu_*(k_^2+((m_^2+1)./(r_.*r_))),v_);
        RHS_box=v_;
    end
    
    %w
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %rho
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %P
    if EigFun_==5
        P_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=-fT_((1i*m_./r_),P_);
    end
    
    %Do not evaluate on the boundaries:
    if and(EigFun_==0,or(T_==1,T_==-1))
        Info_=1;
        return
    end
    %Trim boundaries:
    %We don't wish to evaluate this equation on the boundaries because of
    %the boundary conditions.
    if not(length(LHS_box)==1)
        LHS_box=LHS_box(2:Z_-1,:);
    end
    if not(length(RHS_box)==1)
        RHS_box=RHS_box(2:Z_-1,:);
    end
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NS3:
%Vertical Navier Stokes Equation:
if Eqn_==3
    %u
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %v
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %w
    if EigFun_==3
        w_=CMatrix(Z_,box_width(EigFun_));
        Dw_=w_*dcdx_*DMatrix(box_width(EigFun_));
        DDw_=Dw_*dcdx_*DMatrix(box_width(EigFun_));
        
        LHS_box=nu_*DDw_+fT_((nu_./r_),Dw_)-fT_(1i*m_*Omega_+nu_*(k_^2+((m_^2)./(r_.*r_))),w_);
        RHS_box=w_;
    end
    
    %rho
    if EigFun_==4
        rho_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=-rho_;
    end
    
    %P
    if EigFun_==5
        P_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=-1i*k_*P_;
    end
    
    %Do not evaluate on the boundaries:
    if and(EigFun_==0,or(T_==1,T_==-1))
        Info_=1;
        return
    end
    
    %Trim boundaries:
    %We don't wish to evaluate this equation on the boundaries because of
    %the boundary conditions.
    if not(length(LHS_box)==1)
        LHS_box=LHS_box(2:Z_-1,:);
    end
    if not(length(RHS_box)==1)
        RHS_box=RHS_box(2:Z_-1,:);
    end
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mass:
%Mass Continuity Equation
if Eqn_==4
    
    %u
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %v
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %w
    if EigFun_==3
        w_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=(N_^2)*w_;
    end
    
    %rho
    if EigFun_==4
        rho_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=-fT_(1i*m_*Omega_,rho_);
        RHS_box=rho_;
    end
    
    %P
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Div:
if Eqn_==5
    
    %u
    if EigFun_==1
        u_=CMatrix(Z_,box_width(EigFun_));
        Du_=u_*dcdx_*DMatrix(box_width(EigFun_));
        
        LHS_box=Du_+fT_((1./r_),u_);
    end
    
    %v
    if EigFun_==2
        v_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=fT_(1i*m_./r_,v_);
    end
    
    %w
    if EigFun_==3
        w_=CMatrix(Z_,box_width(EigFun_));
        
        LHS_box=1i*k_*w_;
    end
    
    %rho
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    
    %P
    %Nothing
    %This term does not appear in this equation, so we can leave LHS_box
    %and RHS_box empty.
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Boundary Conditions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%% No-Penetration
%u=0 Inner:
if Eqn_==-1
    %u
    if EigFun_==1
        u_=CMatrix(Z_,box_width(EigFun_));
        LHS_box=u_(1,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%u=0 Outer:
if Eqn_==-2
    %u
    if EigFun_==1
        u_=CMatrix(Z_,box_width(EigFun_));
        LHS_box=u_(Z_,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% No-Slip
%v=0 Inner:
if Eqn_==-3
    %v
    if EigFun_==2
        v_=CMatrix(Z_,box_width(EigFun_));
        LHS_box=v_(1,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v=0 Outer:
if Eqn_==-4
    %v
    if EigFun_==2
        v_=CMatrix(Z_,box_width(EigFun_));
        LHS_box=v_(Z_,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%w=0 Inner:
if Eqn_==-5
    %w
    if EigFun_==3
        w_=CMatrix(Z_,box_width(EigFun_));
        LHS_box=w_(1,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%w=0 Outer:
if Eqn_==-6
    %w
    if EigFun_==3
        w_=CMatrix(Z_,box_width(EigFun_));
        LHS_box=w_(Z_,:);
    end
end

end

function x_ = fT_(a_,b_)
%Combines a Nx1 vector (a_) by a NxM matrix (b_) such that each term of the
%vector multiplies a row of the matrix.
%For example:
%
%      (1)       ( 1, 2, 3, 4)               ( 1, 2, 3, 4)
%   a_=(2).   b_=( 5, 6, 7, 8).   fT_(a_,b_)=(10,12,14,16).
%      (3)       ( 9,10,11,12)               (27,30,33,36)
%
%This process is necessary to ensure that each row of LHS_box or RHS_box
%corresponds to an evaluation of terms at a single co-location point.

x_=bsxfun(@times,a_,b_);

end

function [ C_ ] = CMatrix( Z_height_, Z_width_ )
%CMATRIX Creates a Co-location Matrix.
%Given a vector that lists a truncated Chebyshev expansion of a term F, the
%Co-location matrix will return the value of F at each co-location point.
%
%Z_height_ corresponds to the total number of co-location points (including
%the point corresponding to a_=0).
%
%Z_width_ corresponds to the number of terms in the Chebyshev expansion of
%term F.

T_height_ = Z_height_-1;

C_ = zeros(Z_height_,Z_width_);

for a_ = 1:(Z_height_)    %Row
    for b_ = 1:Z_width_     %Column
        
        j_ = a_-1;  %Row
        k_ = b_-1;  %Column
        
        C_(a_,b_) = cos(k_*j_*pi/T_height_);
    end
end

end

function [ d_ ] = DMatrix( Z_ )
%DMATRIX Writes out the differentiation matrix for a truncated Chebyshev
%expansion.
%Given a vector that lists a truncated Chebyshev expansion of a term F, the
%differentiation matrix will return a truncated Chebyshev expansion of
%dF/dc, where c refers to the parameter space of the co-location points
%(i.e. -1<c<+1).
%To get the truncated Chebyshev expansion corresponding to dF/dx, the
%series must then be multiplied by the constant "dcdx_" (defined in the
%main program body).

d_=zeros(Z_,Z_);

for i_=1:Z_
    for j_=(i_+1):Z_
        if (mod(i_+j_,2))
            if i_==1
                d_(i_,j_)=j_-1;
            else
                d_(i_,j_)=2*(j_-1);
            end
        end
    end
end

end