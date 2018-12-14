%edited 12-8: 
tic
saveit = 1;

diffG = 1;
diffC = 1;                  
simR = 0;

dim = 65;
load(['setsofuncorrelatedlandscapes' num2str(dim)],'all_landscapes')
alpha_weak = .1;
alpha_strong = .9;
Ma = 15;
Mp = 5;
m = 4; %perennial mortality, 1/m*no.inds individuals die every year.
alpha = alpha_strong;
sigsq = 1.5;
beta = 1;
initcond = 20; 
generations = 300; 
%choose how many species are initialized on the landscape?
maxsp = 4; 
%this is for choosing the range of abiotic factors over which to partition
STD = 2.5; 
%habitat breadth
%sigsq = .5;
sigsqR = sigsq;
%microsite carrying capacity (individuals)
k = 3;
%proportion of seeds that stay in the microsite
prop = .5; %when k < 2 you have to use prop = .2 to get one individual to stay in the microsite. 
%sps identities
%annuals = [1,3,5,7]; %don't actually need this for the code, just life.
perennials = [2,4]; %REMEMBER TO CHANGE WITH MAXSP.

%determine what size sampling units to use
divisibles = zeros(1,dim);
for K=1:dim;
divisibles(K) = rem(dim,K);
end
scale = find(divisibles == 0);
scale = scale(2:length(scale)-1);

%averages_initialconditions; %executes this m file
load('baselinevals65.mat','val_R','val_G','val_C')
load('originallandscapes65.mat','originallandscapes')

for t = 1:length(all_landscapes);

landsc = all_landscapes{t,1};
A = landsc{1,1};
B = landsc{1,2};

if simR == 1
R = A;
else
    R_a = val_R(t,1)*ones(dim,dim);
    R_p = val_R(t,2)*ones(dim,dim);
end

if diffG == 1
G = B;
else
    G = val_G(t,1)*ones(dim,dim);
end

if diffC == 1
C = B;  
else
    C = val_C(t,1)*ones(dim,dim);
end


%define types of niche differentiation:
npsim = repelem(0,maxsp); 

%Rrange = ((max(max(R)))-(min(min(R))));
%insetR = Rrange/maxsp;
%rangefordiff = (max(max(R))-insetR)-(min(min(R))+insetR); 
%goby = rangefordiff/(maxsp-1);
%npdiffR = (min(min(R))+insetR):goby:(max(max(R))-insetR);

Crange = ((max(max(C)))-(min(min(C))));
insetC = Crange/maxsp;
rangefordiff = (max(max(C))-insetC)-(min(min(C))+insetC); 
goby = rangefordiff/(maxsp-1);
npdiffC = (min(min(C))+insetC):goby:(max(max(C))-insetC);    

Grange = ((max(max(G)))-(min(min(G))));
insetG = Grange/maxsp;
rangefordiff = (max(max(G))-insetG)-(min(min(G))+insetG); 
goby = rangefordiff/(maxsp-1);
npdiffG = (min(min(G))+insetG):goby:(max(max(G))-insetG);    

%how do species partition resources? 


npR = npsim;

if diffG == 1
npG = npdiffG;
else
    npG = npsim;
end

if diffC == 1
npC = npdiffC;  
else
    npC = npsim;
end


%run the simulation
[endallinds,endtotalinds,endrichness,endfitness,endcomp_inter,endcomp_intra] = simulationv2018(R_a,R_p,C,G,dim,npR,npC,npG,generations,initcond,scale,maxsp,sigsq,sigsqR,Ma,Mp,alpha,beta,m,k,prop,perennials,originallandscapes);


if saveit == 1
        filename = ['simulationv2018_dim' num2str(dim) 'uncorrelatedsets' num2str(t) 'neutRdiffGdiffC_adjd_strongintra_strongdd.mat'];
        save(filename);
end
end

toc
