function color = ColorSelection2( choice )
% ColorSelection: Returns 1x3 matrix with RGB values for color to plot
% based on the target input
%

switch choice
    case 1
        color = [0 0.5 0];
    case 2
        color = [1 0 1];
    case 3
        color = [0 1 1];
    case 4
        color = [1 0 0];
    case 5
        color = [0 1 0];
    case 6
        color = [0 0 1];
    case 7
        color = [0.5 0.5 0.5];
    otherwise
        color = [0.5 0 0.5];
end

end

