function [outputseedind] = reproduction_dispersal(dim,maxsp,total,spsmatrix,Ma,Mp,R_a,R_p,npR,prop,sigsqR,gen)

seedind = zeros(dim,dim,maxsp);  
%go through each occupied microsite:
occupied = find(total>0); %because you have your check calculating the new total from above.
rng(gen*0.1270);
for o = 1:length(occupied); %do this microsite by microsite
    %get the x,y coordinate for the occupied(o) you're on.
    getdim = size(spsmatrix);
    [x,y]=ind2sub(getdim(1:2),occupied(o));
    
    whichsps = find(spsmatrix(x,y,:)); %which sps are reproducing
            
    %calculate their ouput for VR
    fecundity = zeros(1,length(whichsps)); %an empty matrix for # offspring by sps, NOT ACCOUNTING FOR NO INDS OF EACH SPS
    leave = zeros(1,length(whichsps)); %for the next step keep track of how many offspring are dispersing
    for p = 1:length(whichsps); %do for each species in microsite
        %make an if loop here for annuals v. perennials where the only
        %difference is Ma or Mp
        if mod(whichsps(p),2) == 1 %if the species is odd, it's an annual
            fecundity(1,p) = round(Ma*exp(-.5*(R_a(x,y)-(npR(whichsps(p))))^2/sigsqR));
            %CENSUS ANNUALS
            %allinds(x,y,whichsps(p)) = allinds(x,y,whichsps(p))+spsmatrix(x,y,(whichsps(p)));
        else %otherwise it's a perennial
            fecundity(1,p) = round(Mp*exp(-.5*(R_p(x,y)-(npR(whichsps(p))))^2/sigsqR));
            %CENSUS PERENNIALS
            %allinds(x,y,whichsps(p)) = allinds(x,y,whichsps(p))+spsmatrix(x,y,(whichsps(p))); %# individuals present from previous generations, + # in this generation.  For perennials, we don't care whether they're the same or different. 
        end
        
        
        stay = ceil(prop*fecundity(1,p))*spsmatrix(x,y,(whichsps(p)));        
        leave(1,p) = (fecundity(1,p)*spsmatrix(x,y,(whichsps(p))))-stay; %MAKE SURE THESE SUM TO FECUNDITY!!!
        seedind(x,y,whichsps(p))=stay; 
       
    
%DISPERSAL
    if leave(p)>0 %doesn't specify species
    todisperseto = ceil((9-1)*rand(1,leave(p))); %draws the number of which adjacent cells randomly
    %the potential locations for dispersal need to be modified
    %for corners and edges (so that adding individuals to seedind is
    %skipped if they fall off the edge.
    %corners
    if x == 1 && y == 1; %you're in the upper left, and only cells 3 - 5 are available.
        disperse = todisperseto(3<todisperseto & todisperseto<7);
        elseif x == 1 && y == dim; %you're in the upper right
            disperse = todisperseto(5<todisperseto);
        elseif x == dim && y == 1 %you're in the lower left
            disperse = todisperseto(1<todisperseto & todisperseto<5);
        elseif x == dim && y == dim %you're in the lower right
            disperse = todisperseto(todisperseto == 8 | todisperseto<3);
    
   %edges
        elseif x == 1;
            disperse = todisperseto(3<todisperseto);
        elseif x == dim;
            disperse = todisperseto(todisperseto == 8 | todisperseto<5);
        elseif y == 1;
            disperse = todisperseto(1<todisperseto & todisperseto<7);
        elseif y == dim;
            disperse = todisperseto(5<todisperseto & todisperseto<3);
    else
            disperse = todisperseto;
    end
                
        
    %now that we have modified dispersal corners we can update
    %seedind: 
    for d = 1:length(disperse);
        if disperse(d) == 1;
            seedind(x-1,y-1,whichsps(p))=seedind(x-1,y-1,whichsps(p))+1;
            elseif disperse(d) == 2;
                seedind(x-1,y,whichsps(p))=seedind(x-1,y,whichsps(p))+1;
            elseif disperse(d) == 3;
                seedind(x-1,y+1,whichsps(p))=seedind(x-1,y+1,whichsps(p))+1;
            elseif disperse(d) == 4;
                seedind(x,y+1,whichsps(p))=seedind(x,y+1,whichsps(p))+1;                    
            elseif disperse(d) == 5;
                seedind(x+1,y+1,whichsps(p))=seedind(x+1,y+1,whichsps(p))+1;                    
            elseif disperse(d) == 6;
                seedind(x+1,y,whichsps(p))=seedind(x+1,y,whichsps(p))+1;
            elseif disperse(d) == 7;
                seedind(x+1,y-1,whichsps(p))=seedind(x+1,y-1,whichsps(p))+1;
        else 
                seedind(x,y-1,whichsps(p))=seedind(x,y-1,whichsps(p))+1;
        end
    end        
    end
    end
end
outputseedind = seedind;
%outputallinds = allinds;
end
