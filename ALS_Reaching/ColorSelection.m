function color = ColorSelection( target )
% ColorSelection: Returns 1x3 matrix with RGB values for color to plot
% based on the target input
%

switch target
    case -98
        color = [1 1 0];
    case -86
        color = [1 0 1];
    case -64
        color = [0 1 1];
    case -34
        color = [1 0 0];
    case 34
        color = [0 1 0];
    case 64
        color = [0 0 1];
    case 86
        color = [0.5 0.5 0.5];
    otherwise
        color = [0.5 0 0.5];
end

end

