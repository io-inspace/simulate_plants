for b = 1:20
load(['simulationv2018_dim65uncorrelatedsets',num2str(b),'neutRneutGdiffC_adjd_strongintra_strongdd.mat'])
abundgen = sum(endtotalinds,3)/initcond;
figure(b)
clf
plot(abundgen)
legend('a1','p1','a2','p2')
end

for b = 1:20
load(['simulationv2018_dim65uncorrelatedsets',num2str(b),'neutRneutGdiffC_adjd_strongintra_strongdd.mat'])
figure(b)
clf
imagesc(endallinds{1,1}(:,:,1))
colorbar
end

for b = 1:20
load(['simulationv2018_dim65uncorrelatedsets',num2str(b),'neutRneutGdiffC_adjd_strongintra_strongdd.mat'])
figure(b)
clf
imagesc(sum(endallinds{1,1},3))
colorbar
end