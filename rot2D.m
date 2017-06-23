function [xr, yr] =  rot2d(x, y, ang)
    %rotate vectors by geometric angle

   % This routine is part of Rob Hetland's OCTANT package:
    %https://github.com/hetland/octant

    
xr = x.*cos(ang) - y.*sin(ang);
yr = x.*sin(ang) + y.*cos(ang);

