%% Alican Bayındır 26.10.2020
% Project-1 MATH214 
% This script contains the Bisection method, the Newton's method 
% and the Secant method
clear all
clc
% 1st Bisection method

E = @(x) (1/4*pi*(1/36*pi)*10^-9)*((13*(x-(-7)/abs(x-(-7))^3)) ...
    + (9*(x-(-4)/abs(x-(-4))^3)) + (6*(x-11/abs(x-11)^3)) ...
    + (3*(x-14/abs(x-14)^3)));

% The closed interval is [-3,10] so a = -3, b = 10;
x_lower = -3;
x_upper = 10;

%tolerance value that is given by project guide
tol = 10^-9;

%c provides us loop
c = 0; 
i = 1;
FA = E(x_lower);

% Bisection algorithm
while c < 1
    
    p = x_lower + (x_upper - x_lower) / 2;
    FP = E(p);
    
    % The part below is to handle the table
    i_for_table = (1:i);
    p_for_table(i) = p;
    e_for_table(i) = FP;
    
    %if the result found
    if (FP == 0 || (x_upper - x_lower) / 2 < tol)
        fprintf('The bisection is %f!\n', p);
        fprintf('The bisection found %d. iteration! \n', i);
        
        format shortEng
        o = 1;
        while o <= i
            if o == 1
            error(1) = p_for_table(1) - E(0);
            else
            error(o) = p_for_table(o) - p_for_table(o - 1);
            end
        o = o + 1;
        end        
        my_table = table(i_for_table', p_for_table', error','VariableNames',{'Iterations' 'Results' 'Errors'})
        return
    end
    %if not found initialize variables and add 1 to iteration.
    i = i + 1;
    if FA*FP > 0
        x_lower = p;
        FA = FP;
    else
        x_upper = p;
    end
end
    cc = table(i_for_table',p_for_table')
%% 2nd Newton-Raphson method
clear all
clc

E = @(x) (1/4*pi*(1/36*pi)*10^-10)*((13*(x-(-7)/abs(x-(-7))^3)) ...
    + (9*(x-(-4)/abs(x-(-4))^3)) + (6*(x-11/abs(x-11)^3)) ...
    + (3*(x-14/abs(x-14)^3)));

tol = 10^-9;
c = 0;
i = 1;
% p0 is the initial approximation value. I am going to choose
% mid point of [-3,10] interval
p0 = 3.5;
% Derivated version of E
E_derivated = (26249563097690031*sign(p0 - 11))/(19342813113834066795298816* ...
    abs(p0 - 11)^4) - (72385158845145237*sign(p0 + 7))/(38685626227668133590597632* ...
    abs(p0 + 7)^4) - (7158971753915463*sign(p0 + 4))/(9671406556917033397649408* ...
    abs(p0 + 4)^4) + (16704267425802747*sign(p0 - 14))/(19342813113834066795298816* ...
    abs(p0 - 14)^4) + 8219560161902939/38685626227668133590597632;

% Newton - Raphson method's algorithm
while c < 1
    p = p0 - E(p0) / E_derivated; 
    
    % The part below is to handle the table
    i_for_table = (1:i);
    p_for_table(i) = p0; 
    e_for_table(i) = p;
    
        if abs(p - p0) < tol
        fprintf('The procedure succeded in iteration %d. and the value is: %f', i, p);
        
        format shortEng
        o = 1;
        while o <= i
            if o == 1
            error(1) = abs(p_for_table(1) - E(0));
            else
            error(o) = abs(p_for_table(o) - p_for_table(o - 1));
            end
        o = o + 1;
        end        
        my_table = table(i_for_table', p_for_table', error','VariableNames',{'Iterations' 'Results' 'Errors'})
        return
    end
    i = i + 1;
    p0 = p;
end
%% 3rd The Secant Method
clear all
clc

E = @(x) (1/4*pi*(1/36*pi)*10^-10)*((13*(x-(-7)/abs(x-(-7))^3)) ...
    + (9*(x-(-4)/abs(x-(-4))^3)) + (6*(x-11/abs(x-11)^3)) ...
    + (3*(x-14/abs(x-14)^3)));

tol = 10^-9;
c = 0;
i = 2;

p0 = -3;
p1 = 10;

q0 = E(p0);
q1 = E(p1);

%The secant method's algorithm
while c < 1
    p = p1 - q1 * (p1-p0) / (q1-q0);
    

    % The part below is to handle the table
    i_for_table = (1:i);
    p_for_table(i) = p;
    
    if abs(p-p1) < tol
        fprintf('The procedure succeded in iteration %d. and the value is: %f', i, p)
        
        format shortEng
        o = 1;
        while o <= i
            if o == 1
            error(1) = p_for_table(1) - E(0);
            else
            error(o) = p_for_table(o) - p_for_table(o - 1);
            end
        o = o + 1;
        end
        my_table = table(i_for_table', p_for_table', error','VariableNames',{'Iterations' 'Results' 'Errors'})
        return
    end
    i = i + 1;
    p0 = p1;
    q0 = q1;
    p1 = p;
    q1 = E(p);
end