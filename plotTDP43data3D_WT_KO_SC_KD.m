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
% WT = wild type
% KO = knock-out
% SC = scramble
% KD = knock-down

load('Data3D.mat');

% Average intensity data
avgIntWTinsideDAPI = mean(WTintensityInsideDAPI,1);
semIntWTinsideDAPI = std(WTintensityInsideDAPI,1)/sqrt(size(WTintensityInsideDAPI,1));
avgIntKOinsideDAPI = mean(KOintensityInsideDAPI,1);
semIntKOinsideDAPI = std(KOintensityInsideDAPI,1)/sqrt(size(KOintensityInsideDAPI,1));
avgIntSCinsideDAPI = mean(ScrambleintensityInsideDAPI,1);
semIntSCinsideDAPI = std(ScrambleintensityInsideDAPI,1)/sqrt(size(ScrambleintensityInsideDAPI,1));
avgIntKDinsideDAPI = mean(KDintensityInsideDAPI,1);
semIntKDinsideDAPI = std(KDintensityInsideDAPI,1)/sqrt(size(KDintensityInsideDAPI,1));

avgIntWToutsideDAPI = mean(WTintensityOutsideDAPI,1);
semIntWToutsideDAPI = std(WTintensityOutsideDAPI,1)/sqrt(size(WTintensityOutsideDAPI,1));
avgIntKOoutsideDAPI = mean(KOintensityOutsideDAPI,1);
semIntKOoutsideDAPI = std(KOintensityOutsideDAPI,1)/sqrt(size(KOintensityOutsideDAPI,1));
avgIntSCoutsideDAPI = mean(ScrambleintensityOutsideDAPI,1);
semIntSCoutsideDAPI = std(ScrambleintensityOutsideDAPI,1)/sqrt(size(ScrambleintensityOutsideDAPI,1));
avgIntKDoutsideDAPI = mean(KDintensityOutsideDAPI,1);
semIntKDoutsideDAPI = std(KDintensityOutsideDAPI,1)/sqrt(size(KDintensityOutsideDAPI,1));

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

SCratioIntoutVsin = ScrambleintensityOutsideDAPI ./ ScrambleintensityInsideDAPI;
SCratioIntoutVsin(isnan(SCratioIntoutVsin)) = 0;
SCratioIntoutVsin(isinf(SCratioIntoutVsin)) = 0;
avgRatioIntSC = mean(SCratioIntoutVsin,1);
semRatioIntSC = std(SCratioIntoutVsin,1)/sqrt(size(SCratioIntoutVsin,1));

KDratioIntoutVsin = KDintensityOutsideDAPI ./ KDintensityInsideDAPI;
KDratioIntoutVsin(isnan(KDratioIntoutVsin)) = 0;
KDratioIntoutVsin(isinf(KDratioIntoutVsin)) = 0;
avgRatioIntKD = mean(KDratioIntoutVsin,1);
semRatioIntKD = std(KDratioIntoutVsin,1)/sqrt(size(KDratioIntoutVsin,1));


% Average volume occupancy data
avgSzWTinoutDAPI = mean(WTsizeInsideOutsideDAPI,1);
semSzWTinoutDAPI = std(WTsizeInsideOutsideDAPI,1)/sqrt(size(WTsizeInsideOutsideDAPI,1));
avgSzKOinoutDAPI = mean(KOsizeInsideOutsideDAPI,1);
semSzKOinoutDAPI = std(KOsizeInsideOutsideDAPI,1)/sqrt(size(KOsizeInsideOutsideDAPI,1));
avgSzSCinoutDAPI = mean(ScramblesizeInsideOutsideDAPI,1);
semSzSCinoutDAPI = std(ScramblesizeInsideOutsideDAPI,1)/sqrt(size(ScramblesizeInsideOutsideDAPI,1));
avgSzKDinoutDAPI = mean(KDsizeInsideOutsideDAPI,1);
semSzKDinoutDAPI = std(KDsizeInsideOutsideDAPI,1)/sqrt(size(KDsizeInsideOutsideDAPI,1));

% Calculate ratio inside/outside DAPI
WTratioSzoutVsin = WTsizeInsideOutsideDAPI(:,2) ./ WTsizeInsideOutsideDAPI(:,1);
WTratioSzoutVsin(isnan(WTratioSzoutVsin)) = 0;
KOratioSzoutVsin = KOsizeInsideOutsideDAPI(:,2) ./ KOsizeInsideOutsideDAPI(:,1);
KOratioSzoutVsin(isnan(KOratioSzoutVsin)) = 0;
SCratioSzoutVsin = ScramblesizeInsideOutsideDAPI(:,2) ./ ScramblesizeInsideOutsideDAPI(:,1);
SCratioSzoutVsin(isnan(SCratioSzoutVsin)) = 0;
KDratioSzoutVsin = KDsizeInsideOutsideDAPI(:,2) ./ KDsizeInsideOutsideDAPI(:,1);
KDratioSzoutVsin(isnan(KDratioSzoutVsin)) = 0;

%% Plot relative intensity by intensity bin
figure('Name', 'TDP43 intensity distribution', 'units', 'normalized', 'position', [0.1 0.2 0.6 0.7]);
subplot(3,3,1:3);  hold on;
errorbar(0:255, avgIntWTinsideDAPI, semIntWTinsideDAPI, 'k');
errorbar(0:255, avgIntKOinsideDAPI, semIntKOinsideDAPI, 'r');
errorbar(0:255, avgIntSCinsideDAPI, semIntSCinsideDAPI, 'm');
errorbar(0:255, avgIntKDinsideDAPI, semIntKDinsideDAPI, 'b');
xlim([0 255]);
xlabel('Intensity value');
ylabel('Relative frequency');
title('Intranuclear TDP43 intensity distribution');
legend(['WT average' char(177) 'SEM, n=' num2str(size(WTintensityInsideDAPI,1))],...
       ['KO average' char(177) 'SEM, n=' num2str(size(KOintensityInsideDAPI,1))],...
       ['SC average' char(177) 'SEM, n=' num2str(size(ScrambleintensityInsideDAPI,1))],...
       ['KD average' char(177) 'SEM, n=' num2str(size(KDintensityInsideDAPI,1))]);

subplot(3,3,4:6);  hold on;
errorbar(0:255, avgIntWToutsideDAPI, semIntWToutsideDAPI, 'k');
errorbar(0:255, avgIntKOoutsideDAPI, semIntKOoutsideDAPI, 'r');
errorbar(0:255, avgIntSCoutsideDAPI, semIntSCoutsideDAPI, 'm');
errorbar(0:255, avgIntKDoutsideDAPI, semIntKDoutsideDAPI, 'b');
xlim([0 255]);
xlabel('Intensity value');
ylabel('Relative frequency');
title('Extranuclear TDP43 intensity distribution');
legend(['WT average' char(177) 'SEM, n=' num2str(size(WTintensityOutsideDAPI,1))],...
       ['KO average' char(177) 'SEM, n=' num2str(size(KOintensityOutsideDAPI,1))],...
       ['SC average' char(177) 'SEM, n=' num2str(size(ScrambleintensityOutsideDAPI,1))],...
       ['KD average' char(177) 'SEM, n=' num2str(size(KDintensityOutsideDAPI,1))]);

subplot(3,3,7);  hold on;
bar(1, mean(mean(WTratioIntoutVsin,2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
bar(2, mean(mean(KOratioIntoutVsin,2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
bar(3, mean(mean(SCratioIntoutVsin,2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 1],'LineWidth',1.5);
bar(4, mean(mean(KDratioIntoutVsin,2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1.5 );
scatter(ones(1, size(WTratioIntoutVsin,1))*1.1, mean(WTratioIntoutVsin,2), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
scatter(ones(1, size(KOratioIntoutVsin,1))*2.1, mean(KOratioIntoutVsin,2), 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
scatter(ones(1, size(SCratioIntoutVsin,1))*3.1, mean(SCratioIntoutVsin,2), 60, 'o', 'MarkerEdgeColor',[1 0.5 1], 'LineWidth',1);
scatter(ones(1, size(KDratioIntoutVsin,1))*4.1, mean(KDratioIntoutVsin,2), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 1],'LineWidth',1);
errorbar(1, mean(mean(WTratioIntoutVsin,2)), std(mean(WTratioIntoutVsin,2))/sqrt(size(WTratioIntoutVsin,1)), 'k', 'LineWidth',1.5);
errorbar(2, mean(mean(KOratioIntoutVsin,2)), std(mean(KOratioIntoutVsin,2))/sqrt(size(KOratioIntoutVsin,1)), 'r', 'LineWidth',1.5);
errorbar(3, mean(mean(SCratioIntoutVsin,2)), std(mean(SCratioIntoutVsin,2))/sqrt(size(SCratioIntoutVsin,1)), 'm', 'LineWidth',1.5);
errorbar(4, mean(mean(KDratioIntoutVsin,2)), std(mean(KDratioIntoutVsin,2))/sqrt(size(KDratioIntoutVsin,1)), 'b', 'LineWidth',1.5);
xlim([0 5]);
%ylim([0 10]);
ylabel('Extra/intra-nuclear Ratio');
title('All TDP43');
xticklabels({' ', 'WT', 'KO', 'SC', 'KD', ' '});
legend('p vs WT',...
        num2str(ranksum( mean(WTratioIntoutVsin,2), mean(KOratioIntoutVsin,2) ), 3),...
        'p vs SC',...
        num2str(ranksum( mean(SCratioIntoutVsin,2), mean(KDratioIntoutVsin,2) ), 3));
legend('boxoff');
legend('Location', 'northwest');

% subplot(3,3,8);  hold on;
% bar(1, mean(mean(WTratioIntoutVsin(:, 1:129),2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
% bar(2, mean(mean(KOratioIntoutVsin(:, 1:129),2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
% bar(3, mean(mean(SCratioIntoutVsin(:, 1:129),2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 1],'LineWidth',1.5);
% bar(4, mean(mean(KDratioIntoutVsin(:, 1:129),2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1.5 );
% scatter(ones(1, size(WTratioIntoutVsin,1))*1.1, mean(WTratioIntoutVsin(:, 1:129),2), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
% scatter(ones(1, size(KOratioIntoutVsin,1))*2.1, mean(KOratioIntoutVsin(:, 1:129),2), 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
% scatter(ones(1, size(SCratioIntoutVsin,1))*3.1, mean(SCratioIntoutVsin(:, 1:129),2), 60, 'o', 'MarkerEdgeColor',[1 0.5 1], 'LineWidth',1);
% scatter(ones(1, size(KDratioIntoutVsin,1))*4.1, mean(KDratioIntoutVsin(:, 1:129),2), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 1],'LineWidth',1);
% errorbar(1, mean(mean(WTratioIntoutVsin(:, 1:129),2)), std(mean(WTratioIntoutVsin(:, 1:129),2))/sqrt(size(WTratioIntoutVsin,1)), 'k', 'LineWidth',1.5);
% errorbar(2, mean(mean(KOratioIntoutVsin(:, 1:129),2)), std(mean(KOratioIntoutVsin(:, 1:129),2))/sqrt(size(KOratioIntoutVsin,1)), 'r', 'LineWidth',1.5);
% errorbar(3, mean(mean(SCratioIntoutVsin(:, 1:129),2)), std(mean(SCratioIntoutVsin(:, 1:129),2))/sqrt(size(SCratioIntoutVsin,1)), 'm', 'LineWidth',1.5);
% errorbar(4, mean(mean(KDratioIntoutVsin(:, 1:129),2)), std(mean(KDratioIntoutVsin(:, 1:129),2))/sqrt(size(KDratioIntoutVsin,1)), 'b', 'LineWidth',1.5);
% xlim([0 5]);
% %ylim([0 10]);
% ylabel('Extra/intra-nuclear Ratio');
% title('Dim TDP43 (int<=128)');
% xticklabels({' ',  'WT', 'KO', 'SC', 'KD', ' '});
% legend('p vs WT',...
%         num2str(ranksum( mean(WTratioIntoutVsin(:, 1:129),2), mean(KOratioIntoutVsin(:, 1:129),2)), 3),...
%         'p vs SC',...
%         num2str(ranksum( mean(SCratioIntoutVsin(:, 1:129),2), mean(KDratioIntoutVsin(:, 1:129),2)), 3));
% legend('boxoff');
% legend('Location', 'northwest');
    
% subplot(3,3,9);  hold on;
% bar(1, mean(mean(WTintensityInsideDAPI,2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
% bar(2, mean(mean(KOintensityInsideDAPI,2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
% bar(3, mean(mean(ScrambleintensityInsideDAPI,2)), 'FaceColor',[1 1 1],'EdgeColor',[1 0 1],'LineWidth',1.5);
% bar(4, mean(mean(KDintensityInsideDAPI,2)), 'FaceColor',[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1.5 );
%scatter(ones(1, size(WTintensityInsideDAPI,1))*1.1, mean(WTintensityInsideDAPI,1), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
%scatter(ones(1, size(KOintensityInsideDAPI,1))*2.1, mean(KOintensityInsideDAPI,1), 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
%scatter(ones(1, size(SCrambleintensityInsideDAPI,1))*3.1, mean(SCrambleintensityInsideDAPI,1), 60, 'o', 'MarkerEdgeColor',[1 0.5 1], 'LineWidth',1);
%scatter(ones(1, size(KDintensityInsideDAPI,1))*4.1, mean(KDintensityInsideDAPI,1), 60, 'o', 'MarkerEdgeColor',[0.5 0.5 1],'LineWidth',1);
%errorbar(1, mean(WTintensityInsideDAPI,1), std(WTintensityInsideDAPI,1)/sqrt(size(WTintensityInsideDAPI,1)), 'k', 'LineWidth',1.5);
%errorbar(2, mean(KOintensityInsideDAPI,1), std(KOintensityInsideDAPI,1)/sqrt(size(KOintensityInsideDAPI,1)), 'r', 'LineWidth',1.5);
%errorbar(3, mean(ScrambleintensityInsideDAPI,1), std(ScrambleintensityInsideDAPI,1)/sqrt(size(ScrambleintensityInsideDAPI,1)), 'm', 'LineWidth',1.5);
%errorbar(4, mean(KDintensityInsideDAPI,1), std(KDintensityInsideDAPI,1)/sqrt(size(KDintensityInsideDAPI,1)), 'b', 'LineWidth',1.5);
xlim([0 5]);
%ylim([0 2]);
ylabel('Average Intensity');
title('Intranuclear TDP43');
xticklabels({' ', 'WT', 'KO', 'SC', 'KD', ' '});
legend('p vs WT',...
       num2str(ranksum( mean(WTratioIntoutVsin(:, 130:256),2), mean(KOratioIntoutVsin(:, 130:256),2)), 3),...
       'p vs SC',...
       num2str(ranksum( mean(SCratioIntoutVsin(:, 130:256),2), mean(KDratioIntoutVsin(:, 130:256),2)), 3));
legend('boxoff');
legend('Location', 'northwest');

%% Plot volume occupancy
figure('Name', 'TDP43 volume occupancy', 'units', 'normalized', 'position', [0.2 0.2 0.7 0.5]);
subplot(1,2,1); hold on;
scatter(WTsizeInsideOutsideDAPI(:,1), WTsizeInsideOutsideDAPI(:,2), 'ok'); 
scatter(KOsizeInsideOutsideDAPI(:,1), KOsizeInsideOutsideDAPI(:,2), 'or');
scatter(ScramblesizeInsideOutsideDAPI(:,1), ScramblesizeInsideOutsideDAPI(:,2), 'om'); 
scatter(KDsizeInsideOutsideDAPI(:,1), KDsizeInsideOutsideDAPI(:,2), 'ob');
errorbar(avgSzWTinoutDAPI(1), avgSzWTinoutDAPI(2), semSzWTinoutDAPI(2), semSzWTinoutDAPI(2), semSzWTinoutDAPI(1), semSzWTinoutDAPI(1), 'k');
errorbar(avgSzKOinoutDAPI(1), avgSzKOinoutDAPI(2), semSzKOinoutDAPI(2), semSzKOinoutDAPI(2), semSzKOinoutDAPI(1), semSzKOinoutDAPI(1), 'r');
errorbar(avgSzSCinoutDAPI(1), avgSzSCinoutDAPI(2), semSzSCinoutDAPI(2), semSzSCinoutDAPI(2), semSzSCinoutDAPI(1), semSzSCinoutDAPI(1), 'm');
errorbar(avgSzKDinoutDAPI(1), avgSzKDinoutDAPI(2), semSzKDinoutDAPI(2), semSzKDinoutDAPI(2), semSzKDinoutDAPI(1), semSzKDinoutDAPI(1), 'b');

% % Fit with linear functions and plot
% WTfitF = polyfit(WTsizeInsideOutsideDAPI(:,1),WTsizeInsideOutsideDAPI(:,2),1);
% WTfitX = WTsizeInsideOutsideDAPI(:,1);
% WTfitY = WTfitF(1) * WTfitX + WTfitF(2);
% plot(WTfitX,WTfitY,'k-');
% 
% KOfitF = polyfit(KOsizeInsideOutsideDAPI(:,1),KOsizeInsideOutsideDAPI(:,2),1);
% KOfitX = KOsizeInsideOutsideDAPI(:,1);
% KOfitY = KOfitF(1) * KOfitX + KOfitF(2);
% plot(KOfitX, KOfitY,'r-');
% 
% SCfitF = polyfit(ScramblesizeInsideOutsideDAPI(:,1),ScramblesizeInsideOutsideDAPI(:,2),1);
% SCfitX = ScramblesizeInsideOutsideDAPI(:,1);
% SCfitY = SCfitF(1) * SCfitX + SCfitF(2);
% plot(SCfitX, SCfitY,'m-');
% 
% KDfitF = polyfit(KDsizeInsideOutsideDAPI(:,1),KDsizeInsideOutsideDAPI(:,2),1);
% KDfitX = KDsizeInsideOutsideDAPI(:,1);
% KDfitY = KDfitF(1) * KDfitX + KDfitF(2);
% plot(KDfitX, KDfitY,'b-');

xlabel('Intranuclear Volume occupancy (%)');
ylabel('Extranuclear Volume occupancy (%)');
title('TDP43 volume occupancy by compartment');
legend('WT','KO','SC','KD',...
       ['avg' char(177) 'SEM'],...
       ['avg' char(177) 'SEM'],...
       ['avg' char(177) 'SEM'],...
       ['avg' char(177) 'SEM'],...
       'Location', 'northeast');
legend('NumColumns', 2);

subplot(1,2,2); hold on;
bar(1, mean(WTratioSzoutVsin), 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
bar(2, mean(KOratioSzoutVsin), 'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'LineWidth',1.5 );
bar(3, mean(SCratioSzoutVsin), 'FaceColor',[1 1 1],'EdgeColor',[1 0 1],'LineWidth',1.5);
bar(4, mean(KDratioSzoutVsin), 'FaceColor',[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1.5 );
scatter(ones(size(WTratioSzoutVsin))*1.1, WTratioSzoutVsin, 60, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
scatter(ones(size(KOratioSzoutVsin))*2.1, KOratioSzoutVsin, 60, 'o', 'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',1);
scatter(ones(size(SCratioSzoutVsin))*3.1, SCratioSzoutVsin, 60, 'o', 'MarkerEdgeColor',[1 0.5 1], 'LineWidth',1);
scatter(ones(size(KDratioSzoutVsin))*4.1, KDratioSzoutVsin, 60, 'o', 'MarkerEdgeColor',[0.5 0.5 1],'LineWidth',1);
errorbar(1, mean(WTratioSzoutVsin), std(WTratioSzoutVsin)/sqrt(numel(WTratioSzoutVsin)), 'k', 'LineWidth',1.5);
errorbar(2, mean(KOratioSzoutVsin), std(KOratioSzoutVsin)/sqrt(numel(KOratioSzoutVsin)), 'r', 'LineWidth',1.5);
errorbar(3, mean(SCratioSzoutVsin), std(SCratioSzoutVsin)/sqrt(numel(SCratioSzoutVsin)), 'm', 'LineWidth',1.5);
errorbar(4, mean(KDratioSzoutVsin), std(KDratioSzoutVsin)/sqrt(numel(KDratioSzoutVsin)), 'b', 'LineWidth',1.5);
xlim([0 5]);
%ylim([0 1.5]);
ylabel('TDP43 extra/intra-nuclear volume occupancy ratio');
title('TDP43 extra/intra-nuclear ratio');
xticklabels({' ', 'WT', 'KO', 'SC', 'KD', ' '});
legend('p vs WT',...
       num2str(ranksum(WTratioSzoutVsin, KOratioSzoutVsin),3),...
       'p vs SC',...
       num2str(ranksum(SCratioSzoutVsin, KDratioSzoutVsin),3));
legend('boxoff');
legend('Location', 'northwest');


%% Clear temporary variables
clear avg* sem* ans
