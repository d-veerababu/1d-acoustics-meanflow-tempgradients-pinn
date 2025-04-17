function [Ur,Ui] = particle_velocity(X,X0,UV0,parameters_p,xeta4,xeta5)

% Predict the acoustic pressure from the trained parameters
P = model(parameters_p,X);
Pr = (1-X)*UV0(1)+X*UV0(3)/X0(4)+(X0(4)-X).*X.*P(1,:);
Pi = (1-X)*UV0(2)+X*UV0(4)/X0(4)+(X0(4)-X).*X.*P(2,:);

% Calculate gradients of the acoustic pressure
Prx = dlgradient(sum(Pr,'all'),X,'EnableHigherDerivatives',true);
Pix = dlgradient(sum(Pi,'all'),X,'EnableHigherDerivatives',true);

xeta4_r = real(xeta4);
xeta4_i = imag(xeta4);
xeta5_r = real(xeta5);
xeta5_i = imag(xeta5);

% Calculate the particle velocity
Ur = xeta4_r.*Pr-xeta4_i.*Pi+xeta5_r.*Prx-xeta5_i.*Pix;
Ui = xeta4_i.*Pr+xeta4_r.*Pi+xeta5_i.*Prx+xeta5_r.*Pix;

end