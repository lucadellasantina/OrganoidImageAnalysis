%% OrganoidImageAnalysis - Quantify proteins in intracellular compartments
%  Copyright (C) 2019-2020 Luca Della Santina
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% Calculate averages and SEM
load('Data2D.mat');

% Average intensity data
avgIntWTinsideDAPI = mean(WTintensityInsideDAPI,1);
semIntWTinsideDAPI = std(WTintensityInsideDAPI,1)/sqrt(size(WTintensityInsideDAPI,1));
avgIntKOinsideDAPI = mean(KOintensityInsideDAPI,1);
semIntKOinsideDAPI = std(KOintensityInsideDAPI,1)/sqrt(size(KOintensityInsideDAPI,1));

avgIntWToutsideDAPI = mean(WTintensityOutsideDAPI,1);
semIntWToutsideDAPI = std(WTintensityOutsideDAPI,1)/sqrt(size(WTintensityOutsideDAPI,1));
avgIntKOoutsideDAPI = mean(KOintensityOutsideDAPI,1);
semIntKOoutsideDAPI = std(KOintensityOutsideDAPI,1)/sqrt(size(KOintensityOutsideDAPI,1));

% Calculate ratio inside/outside DAPI
WTratioIntoutVsin = WTintensityOutsideDAPI ./ WTintensityInsideDAPI;
WTratioIntoutVsin(isnan(WTratioIntoutVsin)) = 0;
WTratioIntoutVsin(isinf(WTratioIntoutVsin)) = 0;
avgRatioIntWT = mean(WTratioIntoutVsin,1);
semRatioIntWT = std(WTratioIntoutVsin,1)/sqrt(size(WTratioIntoutVsin,1));

KOratioIntoutVsin = KOintensityOutsideDAPI ./ KOintensityInsideDAPI;
KOratioIntoutVsin(isnan(KOratioIntoutVsin)) = 0;
KOratioIntoutVsin(isinf(KOratioIntoutVsin)) = 0;
avgRatioIntKO = mean(KOratioIntoutVsin,1);
semRatioIntKO = std(KOratioIntoutVsin,1)/sqrt(size(KOratioIntoutVsin,1));

% Average volume occupancy data
avgSzWTinoutDAPI = mean(WTsizeInsideOutsideDAPI,1);
semSzWTinoutDAPI = std(WTsizeInsideOutsideDAPI,1)/sqrt(size(WTsizeInsideOutsideDAPI,1));

avgSzKOinoutDAPI = mean(KOsizeInsideOutsideDAPI,1);
semSzKOinoutDAPI = std(KOsizeInsideOutsideDAPI,1)/sqrt(size(KOsizeInsideOutsideDAPI,1));

% Calculate ratio inside/outside DAPI
KOratioSzoutVsin = KOsizeInsideOutsideDAPI(:,2) ./ KOsizeInsideOutsideDAPI(:,1);
KOratioSzoutVsin(isnan(KOratioSzoutVsin)) = 0;
WTratioSzoutVsin = WTsizeInsideOutsideDAPI(:,2) ./ WTsizeInsideOutsideDAPI(:,1);
WTratioSzoutVsin(isnan(WTratioSzoutVsin)) = 0;

%% Plot relative intensity by intensity bin
figure('Name', 'TDP43 intensity distribution', 'units', 'normalized', 'position', [0.1 0.2 0.6 0.7]);
subplot(3,3,1:3);  hold on;
errorbar(0:255, avgIntWTinsideDAPI, semIntWTinsideDAPI, 'k');
errorbar(0:255, avgIntKOinsideDAPI, semIntKOinsideDAPI, 'r');
xlim([0 255]);
xlabel('Intensity value');
ylabel('Relative frequency');
title('Intranuclear TDP43 intensity distribution');
legend(['WT average' char(177) 'SEM, n=' num2str(size(WTintensityInsideDAPI,1))], ['KO average' char(177) 'SEM, n=' num2str(size(KOintensityInsideDAPI,1))]);

subplot(3,3,4:6);  hold on;
errorbar(0:255, avgIntWToutsideDAPI, semIntWToutsideDAPI, 'k');
errorbar(0:255, avgIntKOoutsideDAPI, semIntKOoutsideDAPI, 'r');
xlim([0 255]);
xlabel('Intensity value');
ylabel('Relative frequency');
title('Extranuclear TDP43 intensity distribution');
legend(['WT average'  char(177) 'SEM, n=' num2str(size(WTintensityOutsideDAPI,1))], ['KO average' char(177) 'SEM, n=' num2str(size(KOintensityOutsideDAPI,1))]);

subplot(3,3,7);  hold on;
bar(1, mean(mean(WTratioIntoutVsin,2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
bar(2, mean(mean(KOratioIntoutVsin,2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
scatter(ones(1, size(WTratioIntoutVsin,1))*1.1, mean(WTratioIntoutVsin,2), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
scatter(ones(1, size(KOratioIntoutVsin,1))*2.1, mean(KOratioIntoutVsin,2), 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
errorbar(1, mean(mean(WTratioIntoutVsin,2)), std(mean(WTratioIntoutVsin,2))/sqrt(size(WTratioIntoutVsin,1)), 'k', 'LineWidth',1.5);
errorbar(2, mean(mean(KOratioIntoutVsin,2)), std(mean(KOratioIntoutVsin,2))/sqrt(size(KOratioIntoutVsin,1)),'r', 'LineWidth',1.5);
xlim([0 3]);
%ylim([0 10]);
ylabel('Extra/intra-nuclear Ratio');
title(['All TDP43, p=' num2str(ranksum( mean(WTratioIntoutVsin,2), mean(KOratioIntoutVsin,2) ), 3)]);
%legend('WT', 'KO');
xticklabels({' ', 'WT', 'KO', ' '});


subplot(3,3,8);  hold on;
bar(1, mean(mean(WTratioIntoutVsin(:, 1:129),2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
bar(2, mean(mean(KOratioIntoutVsin(:, 1:129),2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
scatter(ones(1, size(WTratioIntoutVsin,1))*1.1, mean(WTratioIntoutVsin(:, 1:129),2), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
scatter(ones(1, size(KOratioIntoutVsin,1))*2.1, mean(KOratioIntoutVsin(:, 1:129),2), 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
errorbar(1, mean(mean(WTratioIntoutVsin(:, 1:129),2)), std(mean(WTratioIntoutVsin(:, 1:129),2))/sqrt(size(WTratioIntoutVsin,1)), 'k', 'LineWidth',1.5);
errorbar(2, mean(mean(KOratioIntoutVsin(:, 1:129),2)), std(mean(KOratioIntoutVsin(:, 1:129),2))/sqrt(size(KOratioIntoutVsin,1)),'r', 'LineWidth',1.5);
xlim([0 3]);
%ylim([0 10]);
ylabel('Extra/intra-nuclear Ratio');
title(['Dim TDP43, p=' num2str(ranksum( mean(WTratioIntoutVsin(:, 1:129),2), mean(KOratioIntoutVsin(:, 1:129),2)), 3)]);
%legend('WT', 'KO');
xticklabels({' ',  'WT', 'KO', ' '});

subplot(3,3,9);  hold on;
bar(1, mean(mean(WTratioIntoutVsin(:, 130:256),2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
bar(2, mean(mean(KOratioIntoutVsin(:, 130:256),2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
scatter(ones(1, size(WTratioIntoutVsin,1))*1.1, mean(WTratioIntoutVsin(:, 130:256),2), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
scatter(ones(1, size(KOratioIntoutVsin,1))*2.1, mean(KOratioIntoutVsin(:, 130:256),2), 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
errorbar(1, mean(mean(WTratioIntoutVsin(:, 130:256),2)), std(mean(WTratioIntoutVsin(:, 130:256),2))/sqrt(size(WTratioIntoutVsin,1)), 'k', 'LineWidth',1.5);
errorbar(2, mean(mean(KOratioIntoutVsin(:, 130:256),2)), std(mean(KOratioIntoutVsin(:, 130:256),2))/sqrt(size(KOratioIntoutVsin,1)),'r', 'LineWidth',1.5);
xlim([0 3]);
%ylim([0 2]);
ylabel('Extra/intra-nuclear Ratio');
title(['Bright TDP43, p=' num2str(ranksum( mean(WTratioIntoutVsin(:, 130:256),2), mean(KOratioIntoutVsin(:, 130:256),2)), 3)]);
%legend('WT', 'KO');
xticklabels({' ', 'WT', 'KO', ' '});


%% Plot volume occupancy
figure('Name', 'TDP43 volume occupancy', 'units', 'normalized', 'position', [0.2 0.2 0.7 0.5]);
subplot(1,2,1); hold on;
scatter(WTsizeInsideOutsideDAPI(:,1), WTsizeInsideOutsideDAPI(:,2), 'ok'); 
scatter(KOsizeInsideOutsideDAPI(:,1), KOsizeInsideOutsideDAPI(:,2), 'or');
errorbar(avgSzWTinoutDAPI(1), avgSzWTinoutDAPI(2), semSzWTinoutDAPI(2), semSzWTinoutDAPI(2), semSzWTinoutDAPI(1), semSzWTinoutDAPI(1), 'k');
errorbar(avgSzKOinoutDAPI(1), avgSzKOinoutDAPI(2), semSzKOinoutDAPI(2), semSzKOinoutDAPI(2), semSzKOinoutDAPI(1), semSzKOinoutDAPI(1), 'r');

% % Fit with linear functions
% WTfitF = polyfit(WTsizeInsideOutsideDAPI(:,1),WTsizeInsideOutsideDAPI(:,2),1);
% WTfitX = WTsizeInsideOutsideDAPI(:,1);
% WTfitY = WTfitF(1) * WTfitX + WTfitF(2);
% plot(WTfitX,WTfitY,'k-');
% 
% KOfitF = polyfit(KOsizeInsideOutsideDAPI(:,1),KOsizeInsideOutsideDAPI(:,2),1);
% KOfitX = KOsizeInsideOutsideDAPI(:,1);
% KOfitY = KOfitF(1) * KOfitX + KOfitF(2);
% plot(KOfitX, KOfitY,'r-');

xlabel('Intranuclear Volume occupancy (%)');
ylabel('Extranuclear Volume occupancy (%)');
title('TDP43 volume occupancy by compartment');
legend('WT individual organoid', 'KO individual organoid', ['WT average' char(177) 'SEM'], ['KO average' char(177) 'SEM'], 'Location', 'southeast');

subplot(1,2,2); hold on;
bar(1, mean(WTratioSzoutVsin), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
bar(2, mean(KOratioSzoutVsin), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
scatter(ones(size(WTratioSzoutVsin))*1.1, WTratioSzoutVsin, 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
scatter(ones(size(KOratioSzoutVsin))*2.1, KOratioSzoutVsin, 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
errorbar(1, mean(WTratioSzoutVsin), std(WTratioSzoutVsin)/sqrt(numel(WTratioSzoutVsin)), 'k', 'LineWidth',1.5);
errorbar(2,  mean(KOratioSzoutVsin), std(KOratioSzoutVsin)/sqrt(numel(KOratioSzoutVsin)),'r', 'LineWidth',1.5);
xlim([0 3]);
ylim([0 1.5]);
ylabel('TDP43 extra/intra-nuclear volume occupancy ratio');
title(['TDP43 extra/intra-nuclear ratio, p=' num2str(ranksum(WTratioSzoutVsin, KOratioSzoutVsin),3)]);
%legend('WT', 'KO');
xticklabels({' ', ' ', 'WT', ' ', 'KO', ' '});

%% Clear temporary variables
clear avg* sem* ans
