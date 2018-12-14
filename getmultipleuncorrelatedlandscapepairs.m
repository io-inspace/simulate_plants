%edited 10-18

loop = 20;
all_landscapes = cell(loop,1);

for rep = 1:loop

saveit=1;

dim = 65;

%PARAMETERS
dimension = dim;
iterations = 1;
fractaldims = 2.5;
landscapes = cell(iterations,length(fractaldims));
fractalD = zeros(iterations,length(fractaldims));
STD = 2.5; %this is the sd it will be adjusted to if norm = 1
sa = 1; %add wiggle to all points?
norm = 1;
dodiagnostics = 1;
Dcutoff = .1; %how close do you want output landscape D to be to fractaldims?

stopper1 = 0;
stopper2 = 0;

[outA, outB, outcorrelationmat,outfractalD_B,outfractalD_A] = uncorrelatedlandscapes(stopper1, stopper2, dodiagnostics,Dcutoff,dimension,iterations,fractaldims,fractalD,landscapes,STD,sa,norm,dim);
A = outA;
B = outB;
correlationmat = outcorrelationmat;
fractalD_A = outfractalD_A;
fractalD_B = outfractalD_B;
all_landscapes{rep,1} = {A, B, correlationmat,fractalD_B,fractalD_A};
end

if saveit == 1
filename = ['setsofuncorrelatedlandscapes' num2str(dimension) '.mat'];
save(filename);
end


