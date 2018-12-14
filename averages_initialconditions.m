%get spatial averages for abiotically-mediated competitive effect and
%germination probability. Then these will be the two processes that vary
%across simulations, allowing me to partition the variance in species richness
%/spatial variation in sps fitness attributed to neither, one or the other,
%or both.


[val_R,val_G,val_C] = getbaselinevals(dim,all_landscapes,sigsq,Ma,Mp,beta,alpha,maxsp);
if saveit == 1
filename = ['baselinevals' num2str(dim) '.mat'];
save(filename)
end
%load('baselinevals65.mat','val_R','val_G','val_C')

%replicate across initial conditions (i.e. sps placement):

[originallandscapes] = replandscapes(initcond,dim,maxsp);
if saveit == 1
filename = ['originallandscapes' num2str(dim) '.mat'];
save(filename)
end
%load('originallandscapes65.mat','originallandscapes')
