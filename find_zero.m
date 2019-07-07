function [ x_, x_err ] = find_zero(fun,x_0,options_)
%SECANT Homebrew secant fzero method.
%This root-finding algorithm primarily uses the secant method to find a
%value of x for which the function 'fun' returns zero.
%   However, if the secant method suggests a guess that is outside the
%   range of a known sign-change, the algorithm will instead switch to a
%   bisection approach for one iteration before returning to the secant
%   method.
%
%Input Parameters:
% - fun
%   This is a function handle for the operation that we wish to find the
%   root of - i.e. we want to find fun(x)=0.
%
% - x_0
%   This is an initial guess as to the location of the root. This can
%   either be a single value, or a pair of values between which fun changes
%   sign.
%       -If a single value guess is given, find_zero will first attempt to
%       find a sign-change near to the initial guess.
%       -If a pair of values is supplied, find_zero will check that the
%       function fun has opposite signs for the two terms. It will abort if
%       not.
%
% - options_ is a structure term containing various fields:
%    -options_.Display
%      Whether to print out progress on the root-finding algorithm ('on')
%      or not ('off')
%         Note - the display alignment will fail for iterations above 9999.
%
%    -options_.ConstGrowth
%      Whether the function is believed to change positively (pos) of
%      negatively (neg) when x increases. Can also supply (unknown) if not
%      known.
%
%    -options_.MaxIter
%      Set maximum number of iterations of the root-finding algorithm.
%
%    -options_.TolX
%      Maximum tolerance on the final answer of x compared to the actual
%      root location.
%
%    -options_.TolFun
%      Maximum tolerance on the final answer of fun(x) compared to zero.
%

%Default settings (if options_ is not supplied).
display_=0;         %Don't display progress.
look_='both';       %Look in both directions for a sign-change (if only a single guess given).
grow_='unknown';    %By default we don't know whether the given function constantly increases/constantly decreases with increasing x.
Iter_=0;            %Initialise iteration counter to 0.
Max_Iter=5000;      %Maximum number of iterations.
x_Tol=0;            %Tolerance in x for root finding.
f_Tol=0;            %Tolerance in f for root finding.

if (nargin>2)
    %options_ has been supplied- check settings.
    if strcmp(options_.Display,'on')
        %Do display progress.
        %Note - the display alignment will fail for iterations above 9999.
        display_=1;
    end
    if isfield(options_,'ConstGrowth')
        %Allocate whether the function is believed to always change
        %positively (pos) or negatively (neg) when x increases.
        if (strcmp(options_.ConstGrowth,'neg'))||(strcmp(options_.ConstGrowth,'pos')); grow_=options_.ConstGrowth; end
    end
    if isempty(options_.MaxIter)==0
        %Set maximum number of iterations.
        Max_Iter=options_.MaxIter;
    end
    if isempty(options_.TolX)==0
        %Maximum tolerance on the final answer of x compared to the actual
        %root location.
        x_Tol=options_.TolX;
    end
    if isempty(options_.TolFun)==0
        %Maximum tolerance on the final answer of fun(x) compared to zero.
        f_Tol=options_.TolFun;
    end
end
if display_==1
    fprintf('\n');
end

N_x_0=length(x_0);
if N_x_0>1
    %x_0 has been supplied as a pair of terms. This means we expect the
    %sign of fun(x) to change between the two terms.
    
    if N_x_0>2
        error('find_zero.m has been given an illegal guess. Please supply either a single term, or a pair representing a range in which a sign change occurs.')
    end
    x_pre=x_0(1);
    x_now=x_0(2);
    
    f_pre=fun(x_pre);
    f_now=fun(x_now);
    Iter_=Iter_+2;
    
    if display_==1
        fprintf('Initial interval: \n x_0: [%.4g,%.4g], f_0: [%.4g,%.4g] \n', x_pre,x_now,f_pre,f_now);
    end
    
    if (f_pre*f_now)>0
        error('find_zero.m: Supplied interval end-points are not of opposite sign.');
    end
else
    %The guess was a single value. The algorithm now needs to find another
    %value which returns the opposite sign from fun compared to the initial
    %guess.
    if display_==1
        fprintf('Seeking a sign-change interval based on initial guess.\n');
    end
    
    f_0=fun(x_0); %This is the value of f at x=x_0.
    Iter_=Iter_+1;
    
    %If the function is known to always have a positive or negative
    %gradient, then we can constrain the search for a sign-change to only
    %go in one direction.
    if (strcmp(grow_,'pos')&&(f_0<0))||(strcmp(grow_,'neg')&&(f_0>0))
        look_='up';
    elseif (strcmp(grow_,'pos')&&(f_0>0))||(strcmp(grow_,'neg')&&(f_0<0))
        look_='down';
    elseif f_0==0
        %If this is the case, then the initial guess was a valid root. We
        %should return it and exit.
        x_=x_0;
        x_err=0;
        if display_==1
            fprintf('Initial guess gives zero. Exiting.\n');
        end
        return
    end
    
    if display_==1
        fprintf('\t    x_0 \t\t  f(x_0) \n');
        fprintf('%11.4g \t %11.4g\n\n',x_0,f_0);
    end
    %We searches outwards by factors of root-2.
    hunt_=1;
    a_=0;
    n_=0;
    
    %x_above is the latest value of x that we've been looking at in the
    %direction of increasing x.
    %f_above is the value fun(x_above).
    x_above=x_0;
    f_above=f_0; 
    %x_below is the latest value of x that we've been looking at in the
    %direction of decreasing x.
    %f_below is the value fun(x_below).
    x_below=x_0;
    f_below=f_0;
    %Initialised to the x=x_0 values.
    
    if display_==1
        fprintf('f-Count \t\t     x_a \t\t  f(x_a) \t\t     x_b \t\t  f(x_b) \n');
    end
    while (hunt_==1)
        a_=a_+1;
        if mod(a_,2)>0
            n_=n_+1;
            %Increase the Step Width:
            step_=0.02*abs(x_0)*(sqrt(2)^n_);
            
            if not(strcmp(look_,'down')) 
                %ie, look='up' or look='both'. In either case we do want to
                %check increasing x, so we do that here:
                
                %Remember the previous values of x_above and f_above:
                x_above_pre=x_above;
                f_above_pre=f_above;
                
                %The next value of x to look at:
                x_above=x_0+step_;
                
                %Find the new value:
                f_above=fun(x_above);
                Iter_=Iter_+1;
            else
                x_above_pre=x_0;
                x_above=x_0;
                f_above_pre=f_0;
                f_above=f_0;
            end
            
            %Potential Speedup - this if-statement could be moved inside
            %the previous one:
            if f_above*f_0<0
                %If this is the case, we've found our sign-change.
                %Set hunt_=0 n order to leave the while loop:
                hunt_=0;
                
                %The sign-change is presumably between x_above and
                %x_above_pre.
                x_pre=x_above_pre;
                x_now=x_above;
                f_pre=f_above_pre;
                f_now=f_above;
            end
            
            if display_==1
                formatSpec = '%7g \t %11.4g \t %11.4g \t ';
                fprintf(formatSpec,Iter_,x_above,f_above);
            end
        else
            if not(strcmp(look_,'up'))
                %ie, look='down' or look='both'. In either case we do want
                %to check decreasing x, so we do that here.
                
                %Remember the previous values of x_below and f_below:
                x_below_pre=x_below;
                f_below_pre=f_below;
                
                %The next value of x to look at:
                x_below=x_0-step_;
                
                %Find the new value:
                f_below=fun(x_below);
                Iter_=Iter_+1;
            else
                x_below_pre=x_0;
                x_below=x_0;
                f_below_pre=f_0;
                f_below=f_0;
            end
            
            if f_below*f_0<0
                %If this is the case, we've found our sign-change.
                %Set hunt_=0 n order to leave the while loop:
                hunt_=0;
                
                %The sign-change is presumably between x_below and
                %x_below_pre.
                x_pre=x_below_pre;
                x_now=x_below;
                f_pre=f_below_pre;
                f_now=f_below;
            end
            
            if display_==1
                formatSpec = '%11.4g \t %11.4g \n';
                fprintf(formatSpec,x_below,f_below);
            end
        end
        
        if Iter_>Max_Iter
            %Abort. We've gone past our maximum number of iterations.
            fprintf('\n\n find_zero.m: Could not find a sign-change interval from the initial guess. Aborting.\n');
            x_=x_0;
            x_err=0;
            return
        end
    end
    %If we've made it this far, then we have a sign-change interval to work
    %with.
    if display_==1
        formatSpec = '\nSign change interval found between x = %.4g and x = %.4g. \n';
        fprintf(formatSpec,x_pre,x_now);
    end
end

%f_now and f_pre now have different signs.
%Assigning which is which:
if f_pre>0
    x_pos=x_pre;
    x_neg=x_now;
    f_pos=f_pre;
    f_neg=f_now;
else
    x_pos=x_now;
    x_neg=x_pre;
    f_pos=f_now;
    f_neg=f_pre;
end

%Default exit criteria:
if (x_Tol==0)
    x_Tol=abs(x_now-x_pre)/1000;
end

%Perform the hybrid Secant-bisection method:
loop_=1;
i_=0;
x_now = 0.5*(x_neg+x_pos); %First place to investigate is halfway through the range.

if display_==1
    fprintf('\n f-Count \t\t\t  x_ \t\t  f_(x_) \t Method \t\t   x_low \t   x_high \n\n');
end
method_='Initial';

while loop_==1
    %Start of loop:
    i_=i_+1;
    if i_>1
        f_now=fun(x_now);
        Iter_=Iter_+1;
        %Based on the sign of f_now, the sign-change interval will move.
        %(Since we know that x_now is somewhere within the interval).
        if f_now>0
            x_pos=x_now;
            f_pos=f_now;
        else
            x_neg=x_now;
            f_neg=f_neg;
        end
    end
    
    if display_==1
        formatSpec = '%7g \t %11.4g \t %11.4g \t %s \t %11.4g, %11.4g \n';
        fprintf(formatSpec,Iter_,x_now,f_now,method_,x_neg,x_pos);
    end
    
    %Secant Prediction:
    x_next=x_now-f_now*(x_now-x_pre)/(f_now-f_pre);
    
    %Now you need to check whether x_next is within the interval:
    method_='Secant';
    inside_=0;
    if x_pos>x_neg
        if (x_next>x_neg)&&(x_pos>x_next)
            inside_=1;
        end
    elseif (x_next<x_neg)&&(x_pos<x_next)
        inside_=1;
    end
    
    %If x_next isn't within the interval, we do bisection instead:
    if inside_==0
        method_='Bisection';
        x_next = (x_neg+x_pos)/2;
    end
    
    %Advance the secanting:
    x_pre=x_now;
    f_pre=f_now;
    x_now=x_next;
    
    %Exit?
    if Iter_>Max_Iter
        %Iterated too many times.
        loop_=0;
    elseif max(abs(x_now-x_neg),abs(x_now-x_pos))<(x_Tol)
        %We're within x_tolerance.
        if (f_Tol>0)&&(abs(f_pre)<f_Tol)
            %We're also within f_tolerance.
            %Note:this is based on the most recent function evaluation.
            loop_=0;
        elseif (f_Tol==0)
            %We don't care about f_tolerance.
            loop_=0;
        end
    end
    
    %Close the loop.
end

x_=x_now;
x_err=max(abs(x_now-x_neg),abs(x_now-x_pos));

if display_==1
    formatSpec = '\n\n Solution found: x = %.4g, f(x) = %.4g.\n';
    fprintf(formatSpec,x_now,f_now);
    fprintf(' (Displaying 4 significant figures. Error in x is (+/-) %.4g.) \n\n',x_err);
end

end

