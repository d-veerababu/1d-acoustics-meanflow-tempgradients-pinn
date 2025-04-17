function [loss_p,gradients_p] = modelLoss_p(parameters_p,X,X0,UV0,xeta1,xeta2,xeta3)

Q = xeta2./xeta1;
S = xeta3./xeta1;
Q_R = real(Q);
Q_I = imag(Q);
S_R = real(S);
S_I = imag(S);

% Make predictions with the initial conditions.
P = model(parameters_p,X);
U = P(1,:);
V = P(2,:);

% Trial neural network
UG = (1-X)*UV0(1)+X*UV0(3)/X0(4)+(X0(4)-X).*X.*U;
VG = (1-X)*UV0(2)+X*UV0(4)/X0(4)+(X0(4)-X).*X.*V;

% Calculate derivatives with respect to X.
UGx = dlgradient(sum(UG,'all'),X,'EnableHigherDerivatives',true);
VGx = dlgradient(sum(VG,'all'),X,'EnableHigherDerivatives',true);

% Calculate second-order derivatives with respect to X.
UGxx = dlgradient(sum(UGx,'all'),X,'EnableHigherDerivatives',true);
VGxx = dlgradient(sum(VGx,'all'),X,'EnableHigherDerivatives',true);

% Calculate loss.
f_U = UGxx+(Q_R.*UGx-Q_I.*VGx)+(S_R.*UG-S_I.*VG);
f_V = VGxx+(Q_R.*VGx+Q_I.*UGx)+(S_R.*VG+S_I.*UG);

zeroTarget_U = zeros(size(f_U),"like",f_U);
zeroTarget_V = zeros(size(f_V),"like",f_V);
loss_U = l2loss(f_U, zeroTarget_U);
loss_V = l2loss(f_V, zeroTarget_V);

loss_p = loss_U + loss_V;

% Calculate gradients with respect to the learnable parameters.
gradients_p = dlgradient(loss_p,parameters_p);

end
