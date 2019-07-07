function [ m_c, k_c, Re_c, E_c ] = SRI_Find_mk( eta_, mu_, N_, m_guess, k_guess, Re_guess, T_, m_range, k_range )
%Dr Luke Robins 2019 luke.robins@cantab.net
%SRI_FIND_MK Finds the critical horizontal and vertical wavenumbers with
%the smallest critical Reynolds number. This method uses a grid-based
%search over supplied m_ and k_ ranges.
%
%Method:
%   This algorithm assumes that if the system is unstable for a given
%   Reynolds number, it will remain unstable for all larger Reynolds
%   numbers provided that other parameters remain unchanged.
%       This means that the algorithm will not necessarily be well-behaved
%       in the presence of closed unstable domains.
%
%   This method initially optimises for marginal stability at m_=m_guess
%   and k_=k_guess. The corresponding critical Reynolds number Re_c is
%   noted.
%   The method then cycles through each possible combination of m_ and k_
%   from the supplied m_ and k_ ranges.
%       At each possible combination, the algorithm investigates whether
%   Re_c yields an unstable result. If the combination of m_ and k_ is
%   stable, the combination is discarded and the algorithm moves on to the
%   next combination.
%       However, if a given combination of m_ and k_ is unstable, the
%   algorithm then once again optimises for marginal stability, and updates
%   the values of Re_c, m_c and k_c.
%
%       Once the algorithm has finished cycling through every possible m_
%   and k_ combination, then it is assumed that Re_c will have taken the
%   minimum critical Reynolds number for marginal stability across the
%   entire m-k grid, with corresponding values for the critical azimuthal
%   wavenumber m_c and critical vertical wavenumber k_c.
%
%       This result is further improved by optimising k_c using the
%   SRI_Find_k method.
%
%Input Parameters:
% - [eta_,mu_,N_]
%       These 3 inputs are described in detail in SRI_solver.m.
%   They describe radius ratio, rotation ratio and the buoyancy frequency.
%       For this set of parameters we are searching for the critical
%   azimuthal wavenumber m_c and vertical wavenumber k_c which minimise
%   the critical Reynolds number for marginal stability (defined in
%   SRI_Find_Re.m).
%
% - m_guess
%       This is the initial guess for the critical azimuthal wavenumber.
%   This guess should be given as a single integer value.
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
%   SRI_solver.m will instead estimate the appropriate value of T_ to be
%   used for each evaluation of SRI_solver. This estimate will be based on
%   the magnitude of Re_ used for each evaluation.
%
% - m_range
%       This is the range of m_ values that the method will investigate.
%   This should be supplied in the form of a list of integers. If this is
%   not supplied the method will derive a suitable m_range based on the
%   initial guess.
%
% - k_range
%       This is the range of k_ values that the method will investigate.
%   This should be supplied in the form of a list of real numbers. If this
%   is not supplied the method will derive a suitable k_range based on the
%   initial guess.
%
%Output Parameters:
% - m_c
%       This is azimuthal wavenumber which locally minimises the critical
%   Reynolds number for marginal stability.
%
% - k_c
%       This is vertical wavenumber which locally minimises the critical
%   Reynolds number for marginal stability.
%
% - Re_c
%       This is the critical Reynolds number corresponding to m_c, k_c. It
%   will be the local minimum critical Reynolds number for marginal
%   stability.
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

if nargin<7;
    T_=-1;
    %This value of T_ tells the method to dynamically estimate T_ every
    %time that SRI_solver is called. This can significantly cut-down on
    %run-time (since smaller values of Re_ can often be solved with shorted
    %Chebyshev expansions, which are quicker) but direct control over T_
    %may be desired instead.
    %For a full explanation see SRI_solver.m.
end
%If m_range and/or k_range aren't specified, then derive default ranges:
if nargin<8;
    %m_range:
    if mu_<eta_^2
        %Need to account for the possibility of an m_=0 instability:
        m_low=0;
    else
        %For mu_>=eta_^2 then m_=0 is always stable.
        m_low=1;
    end
    if m_guess<3
        %Default range is up to m_=3.
        m_high=3;
    else
        %Otherwise, range up to m_guess+1.
        m_high=m_guess+1;
    end
    m_range=m_low:m_high;
end
if nargin<9;
    %k_range:
    if k_guess<10
        %Default k_range
        k_range=0:0.1:20;
    else
        %Increase the range if the k_guess is higher:
        k_range=0:0.1:(k_guess+10);
    end
end

%Check the initial guess:
[ Re_c, E_c ] = SRI_Find_Re( eta_, mu_, N_, m_guess, k_guess, Re_guess, T_ );
if (Re_c==0)||(real(-1i*E_c)==-999)
    %Failed to find an unstable mode at the initial guess.
    error('SRI_Find_mk Failed to find an unstable mode near the initial guess.');
end
%The algorithm updates m_c and k_c throughout the grid-search, to
%correspond to the current smallest critical Reynolds number:
m_c=m_guess;
k_c=k_guess;

%Setting up the grid search:
L_m=length(m_range);
L_k=length(k_range);
for a_=1:L_m
    for b_=1:L_k
        m_=m_range(a_);
        k_=k_range(b_);
        
        %First check that the current value of Re_c is unstable here:
        [E_c,~,~,~,~,~]=SRI_solver(eta_,mu_,m_,k_,Re_c,N_,T_);
        
        %If we have instability here, then optimise to the new smallest
        %Re_c:
        if real(E_c)>0
            [ Re_c, ~ ] = SRI_Find_Re( eta_, mu_, N_, m_, k_, Re_c, T_ );
            m_c=m_;
            k_c=k_;
        end
    end
end

%We now have a grid co-ordinate with the smallest Re_c for instability.
%Optimise k_c around this grid co-ordinate by using SRI_Find_k:
[ k_c, Re_c, E_c ] = SRI_Find_k( eta_, mu_, N_, m_c, k_c, Re_c, T_ );

end
