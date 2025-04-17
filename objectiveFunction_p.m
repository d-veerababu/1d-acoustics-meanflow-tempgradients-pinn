function [loss_p,gradientsV_p] = objectiveFunction_p(parametersV_p,X,X0,UV0,parameterNames_p,parameterSizes_p,xeta1,xeta2,xeta3)

% Convert parameters to structure of dlarray objects.
parametersV_p = dlarray(parametersV_p);
parameters_p = parameterVectorToStruct(parametersV_p,parameterNames_p,parameterSizes_p);

% Evaluate model loss and gradients.
[loss_p,gradients_p] = dlfeval(@modelLoss_p,parameters_p,X,X0,UV0,xeta1,xeta2,xeta3);

% Return loss and gradients for fmincon.
gradientsV_p = parameterStructToVector(gradients_p);
gradientsV_p = extractdata(gradientsV_p);
loss_p = extractdata(loss_p);

end
