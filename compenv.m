%edited 11-30

%for the biotic component, I want to keep track of how the competitive
%environment is correlated with the abiotic environment. In particular, I
%want to know whether the competitive effect from conspecifics is related to
%the abiotic environment

%to do that, in the last generation you want to know rj in each microsite
%for each sps. But then...do you take the average across species? NO.
%not sure what to do with it.

function [inter_compeffect,intra_compeffect] = compenv(germind,spsmatrix,dim,maxsp,C,npC,alpha,sigsq)
intra_compeffect = NaN(dim,dim,maxsp);
inter_compeffect = NaN(dim,dim,maxsp);
%competitive_ratio = NaN(dim,dim,maxsp);
totalgerm = sum(germind,3);
occugerm = find(totalgerm>0);   
for o = 1:length(occugerm) %do this microsite by microsite
    getdim = [dim dim];
    [x,y]=ind2sub(getdim,occugerm(o)); %get the x,y coords that have germinated seedlings
    %openspace = k-sum(spsmatrix(x,y,:)); %calculate how many new plants can fit (0-k)
    %if openspace > 0 %if there's any space
        %calculate the number of germinated seedlings and which species
        %they are that could potentially go in that microsite
        whichsps = find(germind(x,y,:));
        %if sum(germind(x,y,:)) > openspace 
        %if there are more open spaces than germinated seedlings, then
        %none of the seedlings feel a competitive effect.
        %let's say the compeff is 0 if there are fewer germinating
        %individuals than open spaces. If there are more, then:
        allnumindsps = zeros(1,maxsp);
                    %allnumindsps = zeros(1,10); %WTF IS THIS 10?!?! >:(
                                            %since adults have a competitive 
                                            %effect, get the TOTAL #
                                            %individuals and what species
                                            %they are in the microsite (so
                                            %old perennials and germinated
                                            %seedlings).
         for inds = 1:maxsp;
           allnumindsps(1,inds) = sum(germind(x,y,inds))+sum(spsmatrix(x,y,inds));
         end
 
         allwhichsps = find(allnumindsps);
                
         Rbysps = zeros(1,length(allwhichsps)); %need the competitive effect of all species,
                                                           %including adults
         for bysps = 1:length(allwhichsps);
            Rbysps(1,bysps) = alpha*exp((-.5)*(C(x,y)-(npC(allwhichsps(bysps))))^2/(sigsq));
         end
                        
         for j = 1:length(whichsps); %calculate probability of dying only for germinating species
            intra_compeffect(x,y,whichsps(j)) = allnumindsps(whichsps(j))*(allnumindsps(whichsps(j)));
            compeffect = zeros(1,length(allwhichsps)); %but there's a compeffect of all species
                            for m = 1:length(allwhichsps);
                                if whichsps(j) == allwhichsps(m)
                                    compeffect(m) = 0;
                                else  
                                    compeffect(m) = Rbysps(m)*alpha*allnumindsps(allwhichsps(m));
                                end
                            end  
            inter_compeffect(x,y,whichsps(j)) = allnumindsps(whichsps(j))*sum(compeffect);
            %competitive_ratio(x,y,whichsps(j)) = inter_compeffect(x,y,whichsps(j))/intra_compeffect(x,y,whichsps(j));
         end
end
%normed_competitive_ratio = (competitive_ratio-min(min(min(competitive_ratio))))/(max(max(max(competitive_ratio)))-min(min(min(competitive_ratio))));

%correlation_ratio = cell(1,length(scale));
%correlation_intra = cell(1,length(scale));
%correlation_inter = cell(1,length(scale));
%Copt_proxy = cell(1,maxsp);
   %Cprox = reshape(Copt_proxy{1,n},size(Copt_proxy{1,n},1)^2,1);
    %compsp = reshape(compeff(:,:,n),size(compeff(:,:,n),1)^2,1);
    %corrtoavg_ratio = cell(1,maxsp);
    %corrtoavg_inter = cell(1,maxsp);
    %corrtoavg_intra = cell(1,maxsp);
        %for z = 1:length(scale)
            %for n = 1:maxsp
            %Copt_proxy{1,n} = abs(C-npC(n));
           % numgrid = dim/scale(z);
            %recordcorr_ratio = zeros(numgrid,numgrid);
            %recordcorr_inter = zeros(numgrid,numgrid);
            %recordcorr_intra = zeros(numgrid,numgrid);
            %indexv = 0;
           % for v = 1:numgrid
                %indexw = 0;
                %for w = 1:numgrid
                    %block_prox = Copt_proxy{1,n}((indexv+1):((indexv)+scale(z)),(indexw+1):(w*scale(z)));
                    %block_intrasp = intra_compeffect((indexv+1):((indexv)+scale(z)),(indexw+1):(w*scale(z)),n);                                    
                    %normed_intrasp = (block_intrasp-min(min(block_intrasp)))/(max(max(block_intrasp))-min(min(block_intrasp)));
                    %block_intersp = inter_compeffect((indexv+1):((indexv)+scale(z)),(indexw+1):(w*scale(z)),n);                                    
                    %normed_intersp = (block_intersp-min(min(block_intersp)))/(max(max(block_intersp))-min(min(block_intersp)));                   
                    %block_ratio = normed_competitive_ratio((indexv+1):((indexv)+scale(z)),(indexw+1):(w*scale(z)),n);
                    %I = ~isnan(block_ratio);
                    %if sum(sum(I))==0
                        %recordcorr_ratio(v,w) = NaN;
                        %recordcorr_inter(v,w) = NaN;
                        %recordcorr_intra(v,w) = NaN;
                    %else
                    %recordcorr_ratio(v,w) = corr(block_prox(I),block_ratio(I));
                    %recordcorr_inter(v,w) = corr(block_prox(I),normed_intersp(I));
                    %recordcorr_intra(v,w) = corr(block_prox(I),normed_intrasp(I));
                   % end
                    %indexw = (w*scale(z)); 
                %end
                %indexv = (v*scale(z));    
           % end
            %corrtoavg_ratio{1,n} = recordcorr_ratio;
            %corrtoavg_inter{1,n} = recordcorr_inter;
            %corrtoavg_intra{1,n} = recordcorr_intra;
            %end
        %correlation = cat(3,corrtoavg{:});
        %correlation_ratio{1,z} = corrtoavg_ratio;
        %correlation_inter{1,z} = corrtoavg_inter;
        %correlation_intra{1,z} = corrtoavg_intra;
        %avgcorr{1,z} = nanmean(correlation,3);
        %end
end
