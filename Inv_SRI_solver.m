function [ E_, u_k, v_k, w_k, rho_k, P_k ] = Inv_SRI_solver( eta_, mu_, m_, k_, N_, T_ )
%Inv_SRI_solver
%Dr Luke Robins 2019 luke.robins@cantab.net
%
%   This MatLab program numerically evaluates the inviscid
%Strato-Rotational Instability (SRI) system for any given set of
%parameters.
%
%   Here the Taylor-Couette system is described in cylindrical polar
%co-ordinates [r,theta,z], with the axis of rotation aligned with the
%z-axis, and gravity pointing in the negative-z direction.
%
%   This code should work for any consistent choice of lengthscale and
%timescale. We have typically used it with the gap-width being our choice
%of lengthscale and the reciprocal of the inner rotation rate being our
%choice of timescale. In this choice of units the buoyancy Frequency N_ is
%equivalent to the reciprocal of the Froude Number at the inner radius.
%
%Input Parameters:
%  - eta_
%    This is the radius ratio of the Taylor-Couette system, defined such
%    that eta_ is equal to the inner radius divided by the outer radius.
%
%  - mu_
%    This is the angular velocity ratio of the Taylor-Couette system,
%    defined such that mu_ is equal to the angular velocity of the outer
%    cylinder divided by the angular velocity of the inner cylinder.
%
%  - m_
%    This is any given horizontal (theta) wavenumber of the system. This
%    term must be a real integer.
%
%  - k_
%    This is any given vertical (z) wavenumber of the system.
%
%  - N_
%    This is the Buoyancy Frequency of the vertically-stratified fluid
%    within the Taylor-Couette system. It can be considered to be a measure
%    of the strength of stratification.
%
%  - T_
%    This is the term in the Chebyshev expansion of the system at which we
%    truncate the expansion. The larger this number, the more accurate the
%    evaluation of the system - but also the longer the program will take
%    to run.
%       Note that, because any Chebyshev expansion includes a zeroth term,
%       there will actually be a total T_+1 terms in the Chebyshev
%       expansion. This is sometimes denoted as Z_=T_+1 in the code.
%
%
%Output Parameters:
%  - E_
%    This is the complex eigenvalue corresponding to the mode with the
%    largest real growth rate. The real growth rate is the real component
%    of E_, while the real frequency is equal to the negative of the
%    imaginary component.
%
%  - u_k, v_k, w_k, rho_k, P_k
%    These five vectors are complex eigenvectors corresponding to the
%    perturbation functions of the mode with the largest real growth rate.
%    They can be considered to be Chebyshev term expansions of their
%    respective complex amplitude functions.
%       The first three vectors correspond to the radial (u_k),
%    angular (v_k) and vertical (w_k) velocity perturbations to the basic
%    state flow.
%       The fourth vector (rho_k) corresponds to the density perturbation
%    to the basic state stratification.
%       The fifth vector (P_k) corresponds to the pressure perturbation to
%    the basic state.
%
%       For simplicity, these five perturbation functions are referred to
%    as Eigen-Functions below.
%
%Required Files:
%  - GenEig.m
%    The function in this file utilises a method to phrase any set of
%    one-dimensional linear differential equations with an eigenvalue as a
%    generalised eigenvalue problem.
%
%  - Inv_SRI_equations.m
%    The function in this file describes the viscous SRI system of
%    differential equations in a format that can GenEig.m can read.

function [LHS_box, RHS_box, Info_] = eqns_( Eqn_e, EigFun_e, N_e )
%The MatLab function defined in Inv_SRI_equations.m needs to be redefined
%to fit into the format required by the GenEig.m method.
[LHS_box, RHS_box, Info_] = Inv_SRI_equations( Eqn_e, EigFun_e, N_e, mu_, eta_, m_, k_, N_);
end

%Use GenEig to construct the appropriate matrices:
[ LHS_, RHS_ ] = GenEig( @eqns_, T_ );

%Find the Eigenvectors (Evs_) and Eigenvalues (Es_)
[ Evs_, Es_ ] = eig( LHS_, RHS_ );
%Es_ is a diagonal square matrix where each Es_(i_,i_) is an eigenvalue and
%everything else is 0.
%Evs_ - each column is an eigenvector corresponding to the eigenvalue in
%the same column of Es_.
%   Note that each eigenvector is in fact five Chebyshev expansions
%   concatenated together, to denote the Eigen-Functions:
%   [u_k,v_k,w_k,rho_k,P_k] (in that order).

%Filter out any spurious Eigen-Function expansions:
[ Safe_Es_, Safe_Evs_ ] = Spurious_Filter( Es_, Evs_, T_, @eqns_ );

%Prepared result in case there are no non-spurious results:
E_=-999; u_k=0; v_k=0; w_k=0; rho_k=0; P_k=0;
    %Note that an eigenvector of E_=-999 is always treated as a spurious
    %result. This approach is also used within the GenEig.m method.

if not(isempty(Safe_Es_))
    %Find the eigenfunction with the largest growth rate:
    [ E_, EigVect_ ] = Largest_E( Safe_Es_, Safe_Evs_ );
    
    %Extract the individual variables for that eigenfunction:
    u_k = Eig_Function_Splitter( EigVect_, 1, T_, @eqns_ );
    v_k = Eig_Function_Splitter( EigVect_, 2, T_, @eqns_ );
    w_k = Eig_Function_Splitter( EigVect_, 3, T_, @eqns_ );
    rho_k = Eig_Function_Splitter( EigVect_, 4, T_, @eqns_ );
    P_k = Eig_Function_Splitter( EigVect_, 5, T_, @eqns_ );
    
end

%Print the result to the command line.
fprintf('\n E_ = %6.4g; \t eta_ = %6.4g; \t mu_ = %6.4g; \t m_ = %6.4g; \t k_ = %6.4g; \t N_ = %6.4g; \t T_ = %7g; \n ', real(E_), eta_, mu_, m_, k_, N_, T_ );

end

function [ Safe_Es_, Safe_Evs_ ] = Spurious_Filter( Es_, Evs_, T_, equations_ )
%SPURIOUS_FILTER Scans through a set of Eigen-Functions and removes the
%spurious results. Since each eigenvector is multiple Chebyshev expansions
%concatenated together, this requires that we first split the eigenvector
%into the respective expansions.
%
%   We define a result as spurious if any of the Chebyshev expansions fails
%to exhibit decay as one moves to higher terms in the expansion.
%   This is checked by taking the absolute values of the all terms in the
%expansion. We then compare the sum of the squares of the absolute values
%of the second half of the expansion to that of the first half of the
%expansion. If the second half divided by the first half is not smaller
%than a pre-defined threshold, the result is considered to be spurious.
%   This check is performed on all of the Chebyshev expansions within a
%given result before that result can be considered "safe".

%Determine our threshold for whether to treat a Chebyshev expansion as
%spurious:
limit_=1e-4;
%A larger number here is more 'forgiving'.

%Find out how many Eigen-Functions there are.
[~,~,EigFun_s] = equations_(0,0,0);
%Now that we know, we can use Eig_Function_Splitter to find individual
%Eigen-Functions.

%Determine the total number of eigenvalues:
L_e=length(Es_(1,:));

%Now to start filtering:
%We need to labelling which results are safe to use: 
%   (1 for safe, 0 for spurious)
safe_=ones(1,L_e);

for fi_=1:L_e
    Current_Ev=Evs_(:,fi_); %Selecting out columns.
    Current_E=Es_(fi_,fi_); %Selecting out eigenvalue.
    
    %Quick sanity check:
    if abs(Current_E)==Inf
        safe_(fi_)=0;
    end
    
    %Now proceed through Current_Ev testing each set of Chebyshev
    %expansions:
    EigFun_=1;
    while ((safe_(fi_)==1)&&(EigFun_<EigFun_s+1))
        EigFun_k=Eig_Function_Splitter(Current_Ev,EigFun_,T_,equations_);
        size_=Spurious_Function(EigFun_k);
        
        if (size_>limit_)
            safe_(fi_)=0;
        end
        
        EigFun_=EigFun_+1;
    end
end

%Now that we've got safe_, we know which results are safe and which aren't.
%Selecting them out the safe Eigenvectors and Eigenvalues:
L_s=sum(safe_);
Safe_Es_=zeros(1,L_s);
Safe_Evs_=zeros(L_e,L_s);
s_=1;
for i_=1:L_e
    if safe_(i_)==1
        Safe_Es_(s_)=Es_(i_,i_);
        Safe_Evs_(:,s_)=Evs_(:,i_);
        s_=s_+1;
    end
end

end

function [ EigFun_k ] = Eig_Function_Splitter( EigVect_, EigFun_, T_, equations_ )
%EIG_FUNCTION_SPLITER Extracts a single Eigen-Function from the list of
%Chebyshev coefficients (i.e. the eigenvector, returning the expansion of
%just that Eigen-Function.
Z_=T_+1;
%Z_ is purely for book-keeping. T_ is the term at which we truncate the
%Chebyshev expansion. However, including the zeroth term of the expansion,
%this means that there are actually Z_=T_+1 terms in an expansion truncated
%at this point.

%Need to find the start and finish of each expansion:
EigFun_start=1;

%We iterate through each Eigen-Function in term from the start of the
%vector. If we want the first Eigen-Function, then we don't iterate at all.
%(i.e the loop doesn't even get entered if EigFun_=1).
for i_=1:(EigFun_-1)
    [~,~,C_mod]=equations_(0,i_,0);
    %C_mod here denotes whether a given Eigen-Function has a modification
    %to the number of terms in its expansion. For example, if C_mod=2 for a
    %given Eigen-Function, then that Eigen-Function has T_+2 terms in its
    %respective Chebyshev expansion.
    %These modifications were sometimes necessary to keep the LHS_ and RHS_
    %matrices square.
    EigFun_start=EigFun_start+Z_+C_mod;
end
[~,~,C_mod]=equations_(0,EigFun_,0);
EigFun_end=EigFun_start+T_+C_mod;

%Extract the Chebyshev expansion of the given Eigen-Function:
EigFun_k=EigVect_(EigFun_start:EigFun_end);

end

function [ size_ratio ] = Spurious_Function( xs_ ) 
%SPURIOUS_FUNCTION Examines a single Chebyshev expansion of an
%Eigen-Function. Returns the ratio between the first half squared and the
%second half squared. Absolutes are taken in order to handle complex
%numbers.

%We want to know what half the length of the expansion is.
L_x_c=floor(length(xs_)/2);

Big_sum=0;   %This will be the sum of the squares of the first half of the expansion.
Small_sum=0; %This will be the sum of the squares of the second half of the expansion.
for j_=1:L_x_c
    Big_sum=Big_sum+(abs(xs_(j_))^2);
    Small_sum=Small_sum+(abs(xs_(L_x_c+j_))^2);
end
size_ratio=Small_sum/Big_sum;

end

function [ E_, EigVect_ ] = Largest_E( Safe_Es_, Safe_Evs_ )
%Selects out the result with the largest real part of the eigenvalue - i.e.
%the largest real growth rate.

L_s=length(Safe_Es_);
E_max=real(Safe_Es_(1));
a_E_max=1; % This is the index of the largest growth rate yet seen.
for i_=2:L_s
    if (real(Safe_Es_(i_))>E_max)&&(real(Safe_Es_(i_))~=Inf)
        E_max=real(Safe_Es_(i_));
        a_E_max=i_;
    end
end

E_=Safe_Es_(a_E_max);
EigVect_=Safe_Evs_(:,a_E_max);

end
