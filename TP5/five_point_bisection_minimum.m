function [argmin, minvalue, tried_x, tried_f] = five_point_bisection_minimum(f, low, high)
%FIVE_POINT_BISECTION_MINIMUM find minimum of a function which decrease
%then increase
%   low, high must be int
%   f can be evaluated at integers
% 
%  author: Leman FENG
%  flm8620@gmail.com
tried_x = [];
tried_f = [];

f_low = f(low);
f_high = f(high);
mid = floor((low+high)/2);
f_mid = f(mid);
tried_x = [tried_x, low, high, mid];
tried_f = [tried_f, f_low, f_high, f_mid];

while true
    if mid == low || mid == high
        if f_low < f_high
            argmin = low;
            minvalue = f_low;
        else
            argmin = high;
            minvalue = f_high;
        end
        return
    end
    tried_x = [tried_x, mid];
    tried_f = [tried_f, f_mid];
    if f_mid > f_low || f_mid > f_high
        error('impossible')
    end
    
    l_mid = floor((low+mid)/2);
    r_mid = floor((high+mid)/2);
    if l_mid == low
        f_l_mid = f_low;
    else
        f_l_mid = f(l_mid);
        tried_x = [tried_x, l_mid];
        tried_f = [tried_f, f_l_mid];
    end
    
    if r_mid == mid
        f_r_mid = f_mid;
    else
        f_r_mid = f(r_mid);
        tried_x = [tried_x, r_mid];
        tried_f = [tried_f, f_r_mid];
    end
    
    if f_l_mid >= f_mid && f_r_mid >= f_mid
        low = l_mid;
        f_low = f_l_mid;
        high = r_mid;
        f_high = f_r_mid;
    elseif f_l_mid < f_mid && f_r_mid < f_mid
        error('impossible')
    elseif f_l_mid < f_mid
        high = mid;
        f_high = f_mid;
        mid = l_mid;
        f_mid = f_l_mid;
    elseif f_r_mid < f_mid
        low = mid;
        f_low = f_mid;
        mid = r_mid;
        f_mid = f_r_mid;
    end
end

end

