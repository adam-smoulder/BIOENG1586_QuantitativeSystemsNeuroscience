function cueNum = CueIdentifier( target )
% ColorSelection: Returns 1x3 matrix with RGB values for color to plot
% based on the target input
%

switch target
    case -98
        cueNum = 1;
    case -86
        cueNum = 2;
    case -64
        cueNum = 3;
    case -34
        cueNum = 4;
    case 34
        cueNum = 5;
    case 64
        cueNum = 6;
    case 86
        cueNum = 7;
    otherwise
        cueNum = 8;
end

end

