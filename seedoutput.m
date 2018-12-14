function [numseeds] = seedoutput(maxsp,dim,individuals_in,Ma,Mp,R_a,R_p,npR,sigsqR)

numseeds = zeros(dim,dim,maxsp);
fitls = zeros(dim,dim,maxsp);
for sp = 1:maxsp
    for x = 1:dim
        for y = 1:dim
            if mod(sp,2) == 1 %if the species is odd, it's an annual
            fecundity = round(Ma*exp(-.5*(R_a(x,y)-(npR(sp)))^2/sigsqR));
            else %otherwise it's a perennial
            fecundity = round(Mp*exp(-.5*(R_p(x,y)-(npR(sp)))^2/sigsqR));
            end
            fitls(x,y,sp) = fecundity;
        end
    end
end

for eachsp = 1:maxsp
    numinds = individuals_in(:,:,eachsp);
    numseeds(:,:,eachsp) = numinds.*fitls(:,:,eachsp);
    %normseeds = (numseeds - min(min(numseeds)))/(max(max(numseeds))-(min(min(numseeds))));
end
end