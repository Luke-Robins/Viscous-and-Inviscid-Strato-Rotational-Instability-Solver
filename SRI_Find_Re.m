function [ Re_c, E_c ] = SRI_Find_Re( eta_, mu_, N_, m_, k_, Re_guess, T_, Tol_Re, assume_pos )
%Dr Luke Robins 2019 luke.robins@cantab.net
%SRI_FIND_RE Uses a root-finding algorithm to find the critical Reynolds
%number for marginal stability (ie Real(E_)=0).
%
%Input Parameters
% - [eta_,mu_,N_,m_,k_]
%       These 5 inputs are described in detail in SRI_solver.m.
%   They describe radius ratio, rotation ratio, buoyancy frequency,
%   azimuthal wavenumber and vertical wavenumber.
%   For this set of parameters we are searching for a Reynolds number Re_c
%   at which the system is marginally stable (ie the real growth rate is
%   equal to zero).
%
% - Re_guess
%       This is an initial guess as to the critical Reynolds number. This
%   guess can be a single value, in which case the method will search for a
%   sign-change in the growth rate near to the initial guess.
%       Alternatively, the guess can be supplied as a pair of Reynolds
%   numbers. In this case each Reynolds number should return an oppositely
%   signed growth rate, such that the method can search for the sign-change
%   within the supplied range.
%
% - T_
%       This is the number of terms in the Chebyshev expansion used by
%   SRI_solver.m. If T_ is left unsupplied, or if T_=-1 is given, then
%   SRI_solver.m will instead estimate the appropriate value of T_ to be
%   used for each evaluation of SRI_solver. This estimate will be based on
%   the magnitude of Re_ used for each evaluation.
%
% - Tol_Re
%       This is the tolerance in the critical Reynolds number for the
%   search - i.e. how close the final result need to be to the value at
%   which the sign-change occurs. If this term is left unsupplied, it will
%   default to 0.01.
%
% - assume_pos
%       This input is a flag that controls whether the root-finding method
%   assumes increasing Re_ will correlate to an increasing growth-rate.
%   This is typically the case with the SRI, but not always, and it may be
%   the case that reducing Re_ is necessary to find a region of
%   instability with a positive growth rate.
%       assume_pos=1 corresponds to assuming that increasing Re_ does
%   correlate to an increasing growth-rate.
%       Any other value for the flag makes no assumption.
%
%Output Parameters:
% - Re_c
%       This is the critical Reynolds number of marginal stability, for
%   which the realy growth rate of the mode is equal to zero. This value of
%   Re_c will be accurate to within +/-Tol_Re.
%       Note that if assume_pos=1, then the method will always return the
%   smallest Reynolds number that was found that still returned a positive
%   growth rate. This value of Re_c will still be accurate to within
%   +/-Tol_Re.
%
% - E_c
%       This is the complex eigenvalue from solver.m for the corresponding
%   inputs and for Re_=Re_c.
%       The real component of this term is the growth rate, and for this
%   method the final real value should be very close to zero. If
%   assume_pos=1, then the final real value will be very small but still
%   positive.
%       The imaginary component of this term is the negative of the
%   frequency of the marginally stable mode.
%

%Default values for empty inputs:
if nargin<7
    T_=-1;
    %This value of T_ tells the method to dynamically estimate T_ every
    %time that SRI_solver is called. This can significantly cut-down on
    %run-time (since smaller values of Re_ can often be solved with shorted
    %Chebyshev expansions, which are quicker) but direct control over T_
    %may be desired instead.
    %For a full explanation see SRI_solver.m.
end
if nargin<8
    Tol_Re=0.01;
end
if nargin<9
    assume_pos=1;
end

%We want to record the smallest value of Re that still returned a positive
%growth rate. If we've assumed that increasing Re will correspondingly
%increase the growth rate, then the final value of Re_pos will what we
%return as Re_c.
Re_pos=Inf;
E_pos=-999;

%Define a function handle that takes a Reynolds number and evaluates
%SRI_solver.m, returning the complex growth rate:
    function E_c = f_Ec_Re(Re_)
        %Sanitise the input value of Re_:
        %(Negative or complex values are not desired.)
        Re_=real(Re_);
        if Re_<0
            Re_=0;
        end
        
        %Evaluate:
        [E_c,~,~,~,~,~] = SRI_solver( eta_, mu_, m_, k_, Re_, N_, T_ );
    end

%Define a further function handle that returns the real component of the
%growth rate:
    function E_r = f_Er_Re(Re_)
        E_c = f_Ec_Re(Re_);
        E_r = real(E_c);
        
        %Track smallest Reynolds number with a positive growth rate.
        if (real(E_c)>0)&&(Re_<Re_pos)
            Re_pos=Re_;
            E_pos=E_c;
        end
    end

%Set options for the search:
options_fzero = optimset('TolX',Tol_Re,'Display','off');
if (assume_pos==1); options_fzero.ConstGrowth='pos'; end

%Run the find_zero.m method:
[Re_c] = find_zero(@f_Er_Re,Re_guess,options_fzero);

if assume_pos==1
    %In this case we return the smallest Reynolds number that returned a
    %positive growth rate:
    Re_c=Re_pos;
    E_c=E_pos;
else
    %In this case we still need to retrieve the eigenvalue of marginal
    %stability:
    E_c = f_Ec_Re(Re_c);
end

end