%% to plot the all individual spectra for each subject 


% borrowed from the function osp_plotOverviewSpec.m of Osprey 2.4.0

% Oeltzschner, G., Zöllner, H. J., Hui, S. C., Mikkelsen, M., Saleh, M. G., 
% Tapper, S., & Edden, R. A. (2020). Osprey: Open-source processing, reconstruction & estimation 
% of magnetic resonance spectroscopy data. Journal of neuroscience methods, 343, 108827.

% Created by M.-Y. WANG 04-09-2023



color = [233/255, 109/255, 99/255; 101/255, 153/255, 255/255; 122/255, 186/255, 122/255];

for subi = 1:3
    
   
    GroupString = ['g_' num2str(subi)];
    data = MRSCont.overview.Osprey.sort_data_voxel_1.(GroupString).metab;

     figure(subi)
    plot(data{1}.ppm,real(data{1}.specs(:,1)) ,'color', [color(subi,:)], 'LineWidth', 1); %data
    hold on;
    for kk = 2 : length(data)
        plot(data{kk}.ppm,real(data{kk}.specs(:,1)),'color', [color(subi,:)], 'LineWidth', 1); %data
    end

    ppmRange = [0.5, 4];
    figTitle = 'Individual Specs';
    xlab='';
    ylab = '';
    % Adapt common style for all axes
    set(gca, 'XDir', 'reverse', 'XLim', [ppmRange(1), ppmRange(end)], 'XMinorTick', 'On');
    ticks = get(gca,'XTick');
    set(gca, 'XTick', unique(round(ticks)));
    set(gca, 'LineWidth', 2, 'TickDir', 'out');
    set(gca, 'FontSize', 22);

    % If no y caption, remove y axis
    if isempty(ylab)
        if ~MRSCont.flags.isGUI
            set(gca, 'YColor', 'w');
            % Black axes, white background
            set(gca, 'XColor', 'k');
            set(gca, 'Color', 'w');
            set(gcf, 'Color', 'w');
            %         title(figTitle, 'Interpreter', 'Tex');
        else
            set(gca, 'YColor', MRSCont.colormap.Background);
            set(gca,'YTickLabel',{});
            set(gca,'YTick',{});
            % Dirtywhite axes, light gray background
            set(gca, 'XColor', MRSCont.colormap.Foreground);
            set(gca, 'Color', MRSCont.colormap.Background);
            set(gcf, 'Color', MRSCont.colormap.Background);
            %         title(figTitle, 'Interpreter', 'Tex', 'Color', MRSCont.colormap.Foreground);
        end
    else
        set(gca, 'YColor', 'k');
    end

    box off;
    % xlabel(xlab, 'FontSize', 28);
    % ylabel(ylab, 'FontSize', 28);
end
