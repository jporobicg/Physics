% USAGE: z=zdepth_roms(h,hc,theta,b,N)
% z z-coordinate center of layer (m)
% h depth (m)
% hc = 20;
% theta = 5;
% b = 0.6;
% N number of layers
function z=sigma2zeta(h,hc,theta,b,N)
% Center of layers
ds = 1 / N;
s = - 1 + ds / 2 : ds : 0 - ds / 2;
A = sinh(theta * s) / sinh(theta);
B = (tanh(theta * (s + .5)) - tanh(theta * .5)) / (2 * tanh(theta * .5));
C = (1 - b) * A + b * B;
z = hc * s + (h - hc) * C;
