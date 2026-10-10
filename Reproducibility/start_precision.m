function [one, restore] = start_precision(ndigits)
%START_PRECISION Select arithmetic and restore the previous digits on function exit.
% The calling driver is a function, so onCleanup also works after an error.
if ndigits == 16
    one = 1;
    restore = [];
    return
end
if exist('sym','file') ~= 2 || exist('vpa','file') ~= 2
    error('MARBF:SymbolicRequired', ['This script requires Symbolic Math Toolbox ', ...
        'for its explicitly requested %d-digit calculation. ', ...
        'The numbered table/figure scripts need only base MATLAB.'], ndigits);
end
old = digits(ndigits);
restore = onCleanup(@() digits(old));
one = vpa('1');
end
