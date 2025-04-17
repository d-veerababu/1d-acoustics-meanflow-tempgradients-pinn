function [xeta1,xeta2,xeta3,xeta4,xeta5] = xetas_dl(X,X0,f,p1,T1,m,rho1,u1,gamma,R)
T = T1+m*(X-X0(1));                   % Temperature at x
dTdx = m;                               % Gradient of temperature at x
d2Tdx2 = 0;                             % Double derivative of temperature
c = sqrt(gamma*R*T);                    % Speed of sound at x

a1 = rho1*u1;                           % Coefficient of u^2
a2 = -(p1+rho1*u1^2);                   % Coefficient of u
a3 = p1*u1/T1;                          % Coefficient of u^0

u = (-a2-sqrt(a2^2-4*a1*a3*T))/(2*a1);          % Smallest root of a1*u^2+a2*u1+a3*T = 0                            
u = dlarray(u,"CB");

M = u./c;                                % Mach number at x
dMdx = M.*(1+gamma*M.^2)*dTdx./(2*(1-gamma*M.^2).*T);
k0 = 2*pi*f./c;
alpha = -(dTdx)./(T.*(1-gamma*M.^2));
dudx = -a3*dTdx./(2*a1*u+a2);
d2udx2 = -(2*a1*dudx.^2+a3*d2Tdx2)./(2*a1*u+a2);

p = p1+rho1*u1*(u1-u);
rho = p./(R*T);
dpdx = -a1*dudx;
d2pdx2 = -a1*d2udx2;
beta = (d2pdx2./p)-(2*dpdx*dTdx./(p.*T))+(2*(dTdx./T).^2)-(d2Tdx2./T);

xeta1 = 1-(M.^2)+(2j*M.^2.*dMdx./k0);
xeta2 = -((1-(3+gamma)*M.^2).*alpha+(2j*M.*k0)+(1j*M.*beta./k0)-(2j*M.*alpha.^2./k0));
xeta3 = (k0.^2)+(1j*(2+gamma)*M.*k0.*alpha)-(2j*gamma*k0.*(M.^2).*dMdx)+((2-gamma)*(M.^2).*beta)+((4*gamma-5)*(M.^2).*alpha.^2);

xeta4 = (M./(rho.*c.*k0)).*(k0+1j*M.*alpha*(gamma-1));
xeta5 = (1./(rho.*c.*k0)).*(1j+M.*alpha./k0);
end