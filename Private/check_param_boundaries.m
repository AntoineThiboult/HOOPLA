function check_param_boundaries
%
% Check and plot where calibrated parameters are situated in regards of the
% parameter boundaries. 
%
% The y axis limits represent the low and up parameter boundary. 
% The x axis represents represent the different model parameters. 
% If the calibrated parameters are situtated close to boundary (closer than
% 5% of the paramter range), the point is displayed in red. 
% 
% Can used to check if the paramters boundaries could be possibly be
% narrower/broader
% 
% Please specify the time step below. 

clear
close all
clc

timeStep='3h';

load(fullfile('..','Data',timeStep,'Model_parameters','Calibrated_param.mat'))
nHM=load(fullfile('..','Data',timeStep,'Misc','hydro_model_names.mat'));
nS=load(fullfile('..','Data',timeStep,'Misc','snow_model_names.mat'));
bndHydro=load(fullfile('..','Data',timeStep,'Model_parameters','model_param_boundaries.mat'));
bndSnow=load(fullfile('..','Data',timeStep,'Model_parameters','snow_model_param_boundaries.mat'));

listC=fieldnames(Bestxs);
for iS=1:size(nS.nameS,1)
    figure
    for iC=1:numel(listC)
        for iM=1:size(nHM.nameM,1)
            bx=Bestxs.(listC{iC}).(nHM.nameM{iM,1}); % Best parameter
            lx=cat(2,bndHydro.(nHM.nameM{iM,1}).sMin,bndSnow.(nS.nameS{iS,1}).sMin);% Lower boundary
            ux=cat(2,bndHydro.(nHM.nameM{iM,1}).sMax,bndSnow.(nS.nameS{iS,1}).sMax);% Upper boundary
            
            % Scalling
            sx=(bx-lx)./(ux-lx);
            idsxOut=find(sx>0.95 | sx<0.05); % Close to boundaries points
            sxOut=NaN(1,numel(bx));
            sxOut(idsxOut)=sx(idsxOut);
            
            % Plot
            subplot(4,5,iM)
            plot(sx,'LineStyle','None','Marker','*','Color','b','MarkerSize',4)
            hold on
            plot(sxOut,'LineStyle','None','Marker','*','Color','r','MarkerSize',4) % Highlight in red points close to param bound
            set(gca,'xlim',[0,numel(bx)+1])
            set(gca,'XTick',{},'YTick',{})
            title(nHM.nameM{iM,1})
        end
    end
end