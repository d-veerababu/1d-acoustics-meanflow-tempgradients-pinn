function parameters = buildNet(numLayers,numNeurons)
% Program to build the network
parameters = struct;

sz = [numNeurons 1];
parameters.fc1_Weights = initializeHe(sz,1,"double");
parameters.fc1_Bias = initializeZeros([numNeurons 1],"double");

for layerNumber=2:numLayers-1
    name = "fc"+layerNumber;

    sz = [numNeurons numNeurons];
    numIn = numNeurons;
    parameters.(name + "_Weights") = initializeHe(sz,numIn,"double");
    parameters.(name + "_Bias") = initializeZeros([numNeurons 1],"double");
end

sz = [2 numNeurons];
numIn = numNeurons;
parameters.("fc" + numLayers + "_Weights") = initializeHe(sz,numIn,"double");
parameters.("fc" + numLayers + "_Bias") = initializeZeros([2 1],"double");

end
