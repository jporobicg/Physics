function [u,v] = uv2rho(u,v,angle)

% average vector components to rho points
u = av2(u')';
v = av2(v);

% pad to correct dimension with NaNs at edges
u = u(:,[1 1:end end]);
u(:,[1 end]) = NaN;
v = v([1 1:end end],:);
v([1 end],:) = NaN;

% rotate coordinates
uveitheta = (u+sqrt(-1)*v).*exp(sqrt(-1)*angle);
u = real(uveitheta);
v = imag(uveitheta);
end
