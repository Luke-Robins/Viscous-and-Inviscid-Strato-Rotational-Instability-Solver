function [ LHS_, RHS_ ] = GenEig( equations_, Num_Cheb_terms )
%GENEIG 
%Dr Luke Robins 2019 luke.robins@cantab.net
%
%This file contains a method for approximating a one-dimensional system
%of linear differential equations with a constant unknown growth rate into
%a generalised eigenvalue problem. In this context, the growth rate becomes
%the eigenvalue of the problem.
%
%This requires the equations to be written in a format compatible with the
%method. These are supplied by a separate equations file, denoted by the
%"equations_" function handle in the inputs.
%
%The method allows the degree of accuracy of the approximation to be
%adjusted accordingly. This corresponds to the number of terms in a
%truncated Chebyshev expansion, denoted by the "Num_Cheb_terms" parameter
%in the inputs.
%
%The outputs are a pair of matrices LHS_ and RHS_, corresponding to the
%Left-Hand-Side and Right-Hand-Side of the original equations. Note that
%the equations must be arranged such that the right-hand-side contains all
%terms multiplied by the eigenvalue, and the left-hand-side contains all
%the terms that do not.
%
%       (This method requires that within the linear differential
%       equations, there exist no terms that involve more than one
%       function, nor are there any constant terms that involve no
%       functions. Therefore cycling through each function is sufficient to
%       ensure that all terms in the problem are represented.)
%
%%%%%%%%%%%%
%equations_:
%%%%%%%%%%%%
%   This is a function handle representing the encoded set of generalised
%   eigenvalue equations. It takes 3 inputs:
%       The first input 'Eqn_' is the index of the equation being
%       looked at.
%       The second input 'EigFun_' is the index of the Eigen-Function being
%       looked at.
%       The third input is usually Num_Cheb_terms, the number of terms in
%       the Chebyshev expansion that the solver is using.
%       
%   If the equation Eqn_ doesn't contain a term with an Eigen-Function
%   corresponding to EigFun_, then equations_ returns zero.
%
%   If the equation Eqn_ does contain terms with the Eigen-Function
%   EigFun_, then equations_ returns a pair of 'small' collocation matrices
%   ('LHS_box' and 'RHS_box') representing these terms on the left and
%   right-hand-side of the equation 'Eqn_'.
%       These small collocation matrices can then be assembled by GenEig to
%       produce a larger pair of matrices representing the entire
%       Generalised Eigenvalue Problem in matrix form.
%%%%%%%%%%%
%%%%%%%%%%%
%Need to know information:
%
%equations_ is built to return specific pieces of information regarding the
%generalised eigenvalue problem. These are retrieved by using the following
%non-standard inputs:
%
%  [~,~,Num_Eqns_EFuns] = equations_(0,0,0);
%   This is the number of equations in the system; which is necessarily
%   also the number of Eigen-Functions in the system.
%
%  [~,~,Num_BCs] = equations_(-1,0,0);
%   This is the number of boundary conditions within the system.
%
%
%Further queries:
%   Which equations are evaluated on the boundaries?
%       [~,~,Eval_]=equations(eqn_,0,b_);
%       This will return Eval_=1 if you should not evaluate the equation on
%       the boundary. Otherwise it returns Eval_=0.
%       b_=+/-1; -1 being the inner boundary, +1 being the outer.
%
%   Which Eigen-Functions require variations on their number of
%   coefficients?
%       [~,~,C_Mod]=equations(0,EigFun_,0);
%       This will return C_Mod, the modification to the number of Chebyshev
%       expansion coefficients that the Eigen-Function EigFun_ has. This
%       can be positive or negative; the total is then
%       1+Num_Cheb_terms+C_mod.
%       (The +1 handles the zeroth coefficient.)
%           These modifications exist to ensure that the larger matrix can
%           be kept square, despite any imposed boundary conditions adding
%           further lines.
%           These are typically unnecessary, hence C_mod is generally 0.
%
%References:
%   For further reading on the technique of using spectral methods to
%   convert linear differential problems into generalised eigenvalue
%   problems, see:
%    - "Numerical Analysis of Spectral Methods: Theory and Applications
%    [D. Gottlieb and S. A. Orszag 1977]
%
%    - "Chebyshev and Fourier Spectral Methods"
%    [J. P. Boyd 2001]
%
%    - "Spectral Methods: Fundamentals in Single Domains"
%    [C. Canuto et al. 2006]
%

%Retrieve the number of equations, Eigen-Functions and boundary conditions.
[~,~,Num_Eqns_EFuns] = equations_(0,0,0);
[~,~,Num_BCs] = equations_(-1,0,0);

%Matrix size?
C_Mods=0;   %Total Chebyshev expansion modifications across all Eigen-Functions.
for a_=1:Num_Eqns_EFuns
    [~,~,C_Mod]=equations_(0,a_,0);
    C_Mods=C_Mods+C_Mod;
end
%The width of the matrix will be N+1 for each Eigen-Function, with
%potential Chebyshev expansion modifications:
size_=Num_Eqns_EFuns*(Num_Cheb_terms+1)+C_Mods;
%For the sake of keeping the matrix square, this should be equal to the
%"number of boundary conditions plus the total number of times that each
%equation is evaluated".

%Initialise matrices:
LHS_=zeros(size_,size_);
RHS_=zeros(size_,size_);

%Main equations:
Eq_out=0;
for Eq_i=1:Num_Eqns_EFuns
    %Cycling through the equations.
    
    %Adjusting to the appropriate set of rows:
    %(Each equation is represented by a specific set of rows within the
    %generalised eigenvalue matrices that this program constructs.)
    Eq_in=Eq_out+1;
    [~,~,Eval_in]=equations_(Eq_i,0,-1);
    [~,~,Eval_out]=equations_(Eq_i,0,+1);
    Eq_out=Eq_in+Num_Cheb_terms-(Eval_in+Eval_out);
    
    EF_out=0;
    for EF_j=1:Num_Eqns_EFuns
        %Cycling through the Eigen-Functions.
        
        %Adjusting to the appropriate set of columns:
        %(Each Eigen-Function is represented by a specific set of columns
        %within the generalised eigenvalue matrices.) 
        EF_in=EF_out+1;
        [~,~,C_Mod]=equations_(0,EF_j,0);
        EF_out=EF_in+Num_Cheb_terms+C_Mod;
        
        %Having selected an equation and a Eigen-Function, we have selected
        %a specific region within the left- and right-hand-side generalised
        %eigenvalue matrices.
        
        %Finding inputs:
        [LHS_box,RHS_box,~]=equations_(Eq_i,EF_j,Num_Cheb_terms);
        
        %These are the smaller co-location matrices that need to be
        %inserted into the larger matrices.
        
        %Inserting:
        LHS_(Eq_in:Eq_out,EF_in:EF_out) = LHS_box;
        RHS_(Eq_in:Eq_out,EF_in:EF_out) = RHS_box;
        
    end
end

Eq_in=Eq_out;
%Boundary Conditions:
%   Boundary Conditions are also encoded within equations_ - to select
%   them, negative numbers are used. For example, the second boundary
%   condition can be selected by using equations_(-2,~,~).
for Eq_i=1:Num_BCs
    %Adjusting appropriate row:
    %Boundary Conditions are placed on the bottom set of rows of the
    %matrices.
    Eq_in=1+Eq_in;
    
    EF_out=0;
    for EF_j=1:Num_Eqns_EFuns
        
        %Adjusting appropriate columns:
        %Once again, a specific Eigen-Function has a specific set of
        %columns in the generalised eigenvalue matrices.
        EF_in=EF_out+1;
        [~,~,C_Mod]=equations_(0,EF_j,0);
        EF_out=EF_in+Num_Cheb_terms+C_Mod;
        
        %Finding inputs:
        [LHS_line,RHS_line,~]=equations_(-Eq_i,EF_j,Num_Cheb_terms);
        
        %Inserting:
        LHS_(Eq_in,EF_in:EF_out)=LHS_line;
        RHS_(Eq_in,EF_in:EF_out)=RHS_line;
        
    end
end

%Catching the case where RHS_row=0:
%
%The RHS (Right-Hand-Side) matrix consists of all terms in all equations
%that include a single factor of the eigenvalue. (A linear eigenvalue
%problem never has powers of the eigenvalue other than 1 or 0.)
%
%   However, some equations may not contain the eigenvalue at all. Notably
%   the boundary conditions often don't. In these circumstances, the RHS
%   matrix will have a row that is entirely zero. This can lead to zero
%   eigenvalue results, since zero is technically an eigenvalue of f(x)=0.
%
%   However, zero is also an eigenvalue of interest when it arises
%   'naturally'.
%   
%   In order to avoid these zero-eigenvalues, but nonetheless still enforce
%   f(x)=0, we instead rephrase the problem as f(x)=-999*f(x). This either
%   yields solutions which maintain f(x)=0, or solutions with an eigenvalue
%   of -999, which can safely be ignored.
%
%   The following code searches for any rows where this adjustment is
%   necessary:
for a_=1:size_
    empty_=zeros(1,size_);
    row_=RHS_(a_,:);
    if isequal(row_,empty_)
        
        if isequal(LHS_(a_,:),empty_)
            error('GenEig.m: Completely zero row');
            %This should never happen with a correctly designed equations_
            %file.
        else
            RHS_(a_,:)=LHS_(a_,:);
            LHS_(a_,:)=-999*LHS_(a_,:);
        end
    end
end

end
