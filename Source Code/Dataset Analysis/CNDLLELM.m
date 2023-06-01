path = "Data\06-Apr-2023";
filename = "CN_DL_Lelm_analysis.mat";

datamat = open(fullfile(path, filename))

CN_mid = datamat.CN_DL_Lelm_Mid;
CN_max = datamat.CN_DL_Lelm_Max;
DLs = datamat.DL;

%%
tmpinds = find(CN_max > 1e16)
CN_max(tmpinds) = 0;
%%
figure(1)
subplot(121)
imagesc(log10(CN_max(2:end,:)))
yticks(1:size(CN_max,1)-1)
yticklabels(split(num2str(DLs(2:end)*100)))
xlabel("Subarray size")
ylabel("Diagonal load (%)")
title("Maximum condition number in scene")
c = colorbar();
c.Label.String = "Condition number";


% figure(2)
subplot(122)
imagesc(log10(CN_mid(2:end,:)))
yticks(1:size(CN_mid,1))
yticklabels(split(num2str(DLs*100)))
xlabel("Subarray size")
ylabel("Diagonal load (%)")
title("Condition number between scatterers in scene")
c = colorbar();
c.Label.String = "Condition number";

%%
CNDLMID = figure(3);
plot(DLs, log10(CN_mid(:,21)), 'o')
