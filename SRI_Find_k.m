function [ k_c, Re_c, E_c ] = SRI_Find_k( eta_, mu_, N_, m_, k_guess, Re_guess, T_ )
%Dr Luke Robins 2019 luke.robins@cantab.net
%SRI_FIND_K Uses fminsearch and SRI_Find_Re to find the critical vertical
%wavenumber k_c, for given mu_, eta_, N_, m_, T_. The critical vertical
%wavenumber k_c is the value of the vertical wavenumber corresponding
%to the local minimum in critical Reynolds number for marginal stability.
%Also returns the associated critical Reynolds number.
%
%Input Parameters:
% - [eta_,mu_,N_,m_]
%       These 4 inputs are described in detail in SRI_solver.m.
%   They describe radius ratio, rotation ratio, buoyancy frequency,
%   and azimuthal wavenumber.
%       For this set of parameters we are searching for the critical
%   vertical wavenumber k_c which minimises the critical Reynolds number
%   for marginal stability (defined in SRI_Find_Re.m).
%
% - k_guess
%       This is the initial guess for the critical vertical wavenumber.
%   This guess should be given as a single value.
%
% - Re_guess
%       This is the initial guess for the critical Reynolds number. This
%   guess should be given as a single value. It will be updated after each
%   evaluation of SRI_Find_Re.
%
% - T_
%       This is the number of terms in the Chebyshev expansion used by
%   SRI_solver.m. If T_ is left unsupplied, or if T_=-1 is given, then
%   SRI_Find_Re.m will instead estimate the appropriate value of T_ to be
%   used for each evaluation of SRI_solver. This estimate will be based on
%   the magnitude of Re_ used for each evaluation.
%
%Output Parameters:
% - k_c
%       This is vertical wavenumber which locally minimises the critical
%   Reynolds number for marginal stability.
%
% - Re_c
%       This is the critical Reynolds number corresponding to k_c.
%
% - E_c
%       This is the complex eigenvalue from solver.m for the corresponding
%   inputs with k_=k_c and Re_=Re_c.
%       The real component of this term is the growth rate, and for this
%   method the final real value should be very close to zero. If
%   assume_pos=1, then the final real value will be very small but still
%   positive.
%       The imaginary component of this term is the negative of the
%   frequency of the marginally stable mode.
%

if nargin<7
    T_=-1;
    %This value of T_ tells the method to dynamically estimate T_ every
    %time that SRI_solver is called. This can significantly cut-down on
    %run-time (since smaller values of Re_ can often be solved with shorted
    %Chebyshev expansions, which are quicker) but direct control over T_
    %may be desired instead.
    %For a full explanation see SRI_solver.m.
end

%Tolerance of the search for the critical Reynolds number:
Tol_Re=0.01;

%The concept of a critical vertical wavenumber only makes sense in the
%context where the system becomes unstable from increasing Re_ - i.e.
%locally we expect the growth rate to increase with the Reynolds number:
assume_pos=1;

%Define a function handle for evaluating SRI_Find_Re.m:
    function Re_=f_Re_k(k_)
        %Sanitise the k_ input:
        k_=real(k_);
        
        %Find the critical Reynolds number of marginal stability:
        [Re_,~] = SRI_Find_Re( eta_, mu_, N_, m_, k_, Re_guess, T_, Tol_Re, assume_pos );
        
        %Update Re_guess:
        Re_guess=Re_;
        
        %Sanitise the Reynolds number output (just in case).
        Re_=real(Re_);
    end

%Tolerance of the search for the critical vertical wavenumber:
Tol_k=0.001;
%Set options for the search:
options_fminsearch = optimset('TolX',Tol_k,'MaxIter',30,'MaxFunEvals',30);

%Perform the minimisation search:
[k_c,Re_c] = fminsearch(@f_Re_k,k_guess,options_fminsearch);

%Now recover the corresponding complex eigenvalue:
[E_c,~,~,~,~,~]=SRI_solver(eta_,mu_,m_,k_c,Re_c,N_,T_);

end