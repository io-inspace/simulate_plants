
function [val_R,val_G,val_C] = getbaselinevals(dim,all_landscapes,sigsq,Ma,Mp,beta,alpha,maxsp)
%repro_vals = zeros(length(all_landscapes),2);
%germ_vals = zeros(length(all_landscapes),maxsp);
%factorC_vals = zeros(length(all_landscapes),maxsp);
%for t = 1:length(all_landscapes);
    
%landsc = all_landscapes{t,1};
%A = landsc{1,1};
%B = landsc{1,2};

%you don't need these unless you actually run the diffR code.
%Arange = ((max(max(A)))-(min(min(A))));
%insetA = Arange/maxsp;
%rangefordiff = (max(max(A))-insetA)-(min(min(A))+insetA); 
%goby = rangefordiff/(maxsp-1);
%npdiffA = (min(min(A))+insetA):goby:(max(max(A))-insetA);

%Brange = ((max(max(B)))-(min(min(B))));
%insetB = Brange/maxsp;
%rangefordiff = (max(max(B))-insetB)-(min(min(B))+insetB); 
%goby = rangefordiff/(maxsp-1);
%npdiffB = (min(min(B))+insetB):goby:(max(max(B))-insetB);    

%fecundity_a = zeros(dim,dim,maxsp);
%fecundity_p = zeros(dim,dim,maxsp);
%germination = zeros(dim,dim,maxsp);
%comp = zeros(dim,dim,maxsp);
%for n = 1:maxsp
%    for x = 1:dim
%        for y = 1:dim
%            fecundity_a(x,y,n) = round(Ma*exp(-.5*(A(x,y)-(0))^2/sigsq));
%            %fecundity_a(x,y,n) = round(Ma*exp(-.5*(A(x,y)-(npdiffA(n)))^2/sigsq));
%            fecundity_p(x,y,n) = round(Mp*exp(-.5*(A(x,y)-(0))^2/sigsq));
%            %fecundity_p(x,y,n) = round(Mp*exp(-.5*(A(x,y)-(npdiffA(n)))^2/sigsq));
%            germination(x,y,n) = beta*exp((-.5)*(B(x,y)-(npdiffB(n)))^2/(sigsq));
%            comp(x,y,n) = alpha*exp((-.5)*(B(x,y)-(npdiffB(n)))^2/(sigsq));
%        end
%    end
    
%    syms r c
%    eqn = Ma*exp(-.5*(r-(0))^2/sigsq) == mean2(fecundity_a(:,:,n));
%    val_a = double(solve(eqn,r));
%    eqn = Mp*exp(-.5*(r-(0))^2/sigsq) == mean2(fecundity_p(:,:,n));
%    val_p = double(solve(eqn,r));
%    %eqn = beta*exp((-.5)*(c-(npdiffB(n)))^2/(sigsq)) == mean2(germination(:,:,n));
%    %val_g = double(solve(eqn,c));
%    eqn = alpha*exp((-.5)*(c-(npdiffB(n)))^2/(sigsq)) == mean2(comp(:,:,n));
%    val_c = double(solve(eqn,c));
%    repro_vals(t,1) = val_a(1);
%    repro_vals(t,2) = val_p(1);
%    %germ_vals(t,n) = val_g(1);
%    factorC_vals(t,n) = val_c(1);%,val_p(1),val_g(1),val_c(1)];  
    
%end
%end

val_R = zeros(length(all_landscapes),2);
val_G = zeros(length(all_landscapes),1);
val_C = zeros(length(all_landscapes),1);
for t = 1:length(all_landscapes);
    
landsc = all_landscapes{t,1};
A = landsc{1,1};
B = landsc{1,2};

%you don't need these unless you actually run the diffR code.
%Arange = ((max(max(A)))-(min(min(A))));
%insetA = Arange/maxsp;
%rangefordiff = (max(max(A))-insetA)-(min(min(A))+insetA); 
%goby = rangefordiff/(maxsp-1);
%npdiffA = (min(min(A))+insetA):goby:(max(max(A))-insetA);

Brange = ((max(max(B)))-(min(min(B))));
insetB = Brange/maxsp;
rangefordiff = (max(max(B))-insetB)-(min(min(B))+insetB); 
goby = rangefordiff/(maxsp-1);
npdiffB = (min(min(B))+insetB):goby:(max(max(B))-insetB);    

fecundity_a = zeros(dim,dim,maxsp);
fecundity_p = zeros(dim,dim,maxsp);
germination = zeros(dim,dim,maxsp);
comp = zeros(dim,dim,maxsp);
for n = 1:maxsp
    for x = 1:dim
        for y = 1:dim
            fecundity_a(x,y,n) = round(Ma*exp(-.5*(A(x,y)-(0))^2/sigsq));
            %fecundity_a(x,y,n) = round(Ma*exp(-.5*(A(x,y)-(npdiffA(n)))^2/sigsq));
            fecundity_p(x,y,n) = round(Mp*exp(-.5*(A(x,y)-(0))^2/sigsq));
            %fecundity_p(x,y,n) = round(Mp*exp(-.5*(A(x,y)-(npdiffA(n)))^2/sigsq));
            germination(x,y,n) = beta*exp((-.5)*(B(x,y)-(npdiffB(n)))^2/(sigsq));
            comp(x,y,n) = alpha*exp((-.5)*(B(x,y)-(npdiffB(n)))^2/(sigsq));
        end
    end
end
    
    syms r c g
    spsmeans = mean(fecundity_a,3);
    totalmean = mean2(spsmeans);
    eqn = Ma*exp(-.5*(r-(0))^2/sigsq) == totalmean;
    val_a = double(solve(eqn,r));
    
    spsmeans = mean(fecundity_p,3);
    totalmean = mean2(spsmeans);
    eqn = Mp*exp(-.5*(r-(0))^2/sigsq) == totalmean;
    val_p = double(solve(eqn,r));
    val_R(t,1) = val_a(1);
    val_R(t,2) = val_p(1);
    
    spsmeans = mean(germination,3);
    totalmean = mean2(spsmeans);
    eqn = beta*exp((-.5)*(g-(0))^2/(sigsq)) == totalmean;
    val_g = double(solve(eqn,g));
    val_G(t,1) = val_g(1);
    
    
    spsmeans = mean(comp,3);
    totalmean = mean2(spsmeans);
    eqn = alpha*exp((-.5)*(c-(0))^2/(sigsq)) == totalmean;
    val_c = double(solve(eqn,c));
    val_C(t,1) = val_c(1);
   
end



%vals = zeros(length(all_landscapes),2);
%for t = 1:length(all_landscapes);
    
%landsc = all_landscapes{t,1};
%A = landsc{1,1};
%B = landsc{1,2};
%fecundity_a = zeros(dim,dim);
%fecundity_p = zeros(dim,dim);
%germination = zeros(dim,dim);
%comp = zeros(dim,dim);
%    for x = 1:dim
%        for y = 1:dim
%            fecundity_a(x,y) = round(Ma*exp(-.5*(A(x,y)-(0))^2/sigsq));
%            fecundity_p(x,y) = round(Mp*exp(-.5*(A(x,y)-(0))^2/sigsq));
%            germination(x,y) = beta*exp((-.5)*(B(x,y)-(0))^2/(sigsq));
%            comp(x,y) = alpha*exp((-.5)*(B(x,y)-(0))^2/(sigsq));
%        end
%    end
%    syms r c
%    eqn = Ma*exp(-.5*(r-(0))^2/sigsq) == mean2(fecundity_a);
%    val_a = double(solve(eqn,r));
%    eqn = Mp*exp(-.5*(r-(0))^2/sigsq) == mean2(fecundity_p);
%    val_p = double(solve(eqn,r));
%    eqn = beta*exp((-.5)*(c-(0))^2/(sigsq)) == mean2(germination);
%    val_g = double(solve(eqn,c));
%    eqn = alpha*exp((-.5)*(c-(0))^2/(sigsq)) == mean2(comp);
%    val_c = double(solve(eqn,c));
%    vals(t,1) = mean([val_a(1),val_p(1)]);
%    vals(t,2) = mean([val_g(1),val_c(1)]);   
%end
end

