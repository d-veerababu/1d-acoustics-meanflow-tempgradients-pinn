function [p,dpdx] = varying_temp_flow_bvp1(xtest,x1,m,gamma,R,rho1,u1,p1,T1,f,U0)

solinit = bvpinit(xtest,@guess);
sol = bvp4c(@(x,y)bvpfcn(x,y,U0,x1,m,gamma,R,rho1,u1,p1,T1,f), @(ya,yb)bcfcn(ya,yb,U0,x1,m,gamma,R,rho1,u1,p1,T1,f), solinit);
p = sol.y(1,:);
dpdx = sol.y(2,:);
end

function g = guess(x) % initial guess for y and y'
g = [sin(x);cos(x)];
end

function res = bcfcn(ya,yb,U0,x1,m,gamma,R,rho1,u1,p1,T1,f)
res = [ya(1)-U0(1); yb(1)-U0(2)];
end

function dydx = bvpfcn(x,y,U0,x1,m,gamma,R,rho1,u1,p1,T1,f)
T = T1+m*(x-x1);                    % Temperature at x
dTdx = m;                           % Gradient of temperature at x
d2Tdx2 = 0;
c = sqrt(gamma*R*T);                % Speed of sound at x

a1 = rho1*u1;                       % Coefficient of u^2
a2 = -(p1+rho1*u1^2);               % Coefficient of u
a3 = p1*u1/T1;                      % Coefficient of u^0

r = roots([a1 a2 a3*T]);            % Roots of a1*u^2+a2*u1+a3*T = 0
% u = r(r>=0)                        % Positive root (Mean flow velocity at x)
u = min(r);                        % Minimum root (Mean flow velocity at x)

M = u/c;                            % Mach number at x
dMdx = M*(1+gamma*M^2)*dTdx/(2*(1-gamma*M^2)*T);
k0 = 2*pi*f/c;
alpha = -(dTdx)/(T*(1-gamma*M^2));
dudx = -a3*dTdx/(2*a1*u+a2);
d2udx2 = -(2*a1*dudx^2+a3*d2Tdx2)/(2*a1*u+a2);

p = p1+rho1*u1*(u1-u);
% rho = p/(R*T);
dpdx = -a1*dudx;
d2pdx2 = -a1*d2udx2;
beta = (d2pdx2/p)-(2*dpdx*dTdx/(p*T))+(2*(dTdx/T)^2)-(d2Tdx2/T);

xeta1 = 1-(M.^2)+(2j*M.^2.*dMdx./k0);
xeta2 = -((1-(3+gamma)*M.^2).*alpha+(2j*M.*k0)+(1j*M.*beta./k0)-(2j*M.*alpha.^2./k0));
xeta3 = (k0.^2)+(1j*(2+gamma)*M.*k0.*alpha)-(2j*gamma*k0.*(M.^2).*dMdx)+((2-gamma)*(M.^2).*beta)+((4*gamma-5)*(M.^2).*alpha.^2);
% [x rho]
dydx = [y(2); -(xeta2./xeta1).*y(2)-(xeta3./xeta1).*y(1)];
end