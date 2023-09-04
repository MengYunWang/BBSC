%% to plot the meanfit and mean spectra 


% borrowed from the function osp_plotMeansSpec.m of Osprey 2.4.0

% Oeltzschner, G., Zöllner, H. J., Hui, S. C., Mikkelsen, M., Saleh, M. G., 
% Tapper, S., & Edden, R. A. (2020). Osprey: Open-source processing, reconstruction & estimation 
% of magnetic resonance spectroscopy data. Journal of neuroscience methods, 343, 108827.

% Created by M.-Y. WANG 04-09-2023


spec = 'metab';
basis = 1;
ind = 1;
sort_fit = 'sort_models_voxel_1';
ppmRange = [0.5, 4];
shift = 0;
[~,min_ppm_index] = min(ppm);
which_spec = spec;
color = [233/255, 109/255, 99/255; 101/255, 153/255, 255/255; 122/255, 186/255, 122/255];
colorFit  = [51/255, 51/255, 51/255];


figTitle = ['mean data \pm SD & mean model: ' spec];
xlab='';
ylab = '';


for subi = 1

    GroupString = ['g_' num2str(subi)];

    fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_fit_' spec])(basis,:,ind);
    fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_fit_' spec])(basis,:,ind);
    data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec])(basis,:,ind);
    data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec])(basis,:,ind);
    baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec])(basis,:,ind);
    baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec])(basis,:,ind);
    residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec])(basis,:,ind);
    residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec])(basis,:,ind);
    ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec])(basis,:,ind);

    if strcmp(spec,'mm')
        MM_clean_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_MM_clean_' spec])(basis,:,ind);
        MM_clean_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_MM_clean_' spec])(basis,:,ind);
    end

    if MRSCont.opts.fit.fitMM
        MM_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_fittMM_' spec])(basis,:,ind);
        MM_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_fittMM_' spec])(basis,:,ind);
    end

    %Calculate SD shadows
    if length(data_sd) > 1
        data_yu = data_mean + data_sd;
        data_yl = data_mean - data_sd;
    end

    if exist('fit_mean', 'var') && length(fit_sd) > 1
        fit_yu = fit_mean + fit_sd;
        fit_yl = fit_mean - fit_sd;
    end
    if exist('baseline_mean', 'var') && length(baseline_sd) > 1
        baseline_yu = baseline_mean + baseline_sd;
        baseline_yl = baseline_mean - baseline_sd;
    end
    if exist('residual_mean', 'var') && length(residual_sd) > 1
        residual_yu = residual_mean + residual_sd;
        residual_yl = residual_mean - residual_sd;
    end
    if exist('MM_clean_sd', 'var') && length(MM_clean_sd) > 1
        MM_clean_yu = MM_clean_mean + MM_clean_sd;
        MM_clean_yl = MM_clean_mean - MM_clean_sd;
    end
    if exist('MM_mean', 'var') && length(MM_sd) > 1 && ~isnan(MM_mean(1))
        MM_yu = MM_mean + MM_sd;
        MM_yl = MM_mean - MM_sd;
    end
    if exist('data_yu', 'var')
        maxshift = max(data_yu);
        maxshift_abs = max(abs(data_yu));
    else
        maxshift = max(data_mean);
        maxshift_abs = max(abs(data_mean));
    end
    shift = maxshift_abs * shift;


    %%% 3. SET UP FIGURE LAYOUT %%%
    % Generate a new figure and keep the handle memorized
    out = figure('Visible','on');
    hold on

    %%% 4. PLOT DATA, FIT, RESIDUAL, BASELINE %%%
    % Determine a positive stagger to offset data, fit, residual, and
    % baseline from the individual metabolite contributions


    if exist('fit_mean', 'var')

        plot(ppm, (zeros(1,length(ppm)) + max(maxshift)), 'Color','k', 'LineWidth', 2); % Maximum Data
        plot(ppm, (zeros(1,length(ppm)) + max(maxshift + abs(min(residual_mean)))), 'Color','k', 'LineStyle','--', 'LineWidth', 1); % Zeroline Residue
        plot(ppm, (zeros(1,length(ppm)) + max(maxshift + abs(min(residual_mean))) + abs(max(residual_mean))), 'Color','k', 'LineWidth', 2); % Max Residue
%         text(ppm(min_ppm_index)-0.05, 0 + max(maxshift +  abs(min(residual_mean))), '0', 'FontSize', 16,'Color','k'); %Zeroline Residual text
        text(ppm(min_ppm_index)-0.05, (0 +max(maxshift + abs(min(residual_mean))) + abs(max(residual_mean))), [num2str(100/max(fit_mean)*abs(max(residual_mean)),'%10.1f') '%'], 'FontSize', 16,'Color','k'); %Max Residue text

    end

    if exist('residual_mean', 'var')
        if length(residual_sd) > 1
            fill([ppm fliplr(ppm)], [(residual_yu+shift+ max(maxshift + abs(min(residual_mean)))) (fliplr(residual_yl)+shift+ max(maxshift + abs(min(residual_mean))))], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow res
        end
    end

    if exist('data_yu', 'var')
        fill([ppm fliplr(ppm)], [data_yu+shift fliplr(data_yl)+shift], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow data
    end

    plot(ppm,data_mean+shift ,'color',[color(subi,:)], 'LineWidth', 3); % Data

    if exist('fit_mean', 'var')
        plot(ppm,fit_mean+shift ,'color', colorFit, 'LineWidth', 2); %Fit

        plot(ppm,residual_mean+shift+ max(maxshift + abs(min(residual_mean))) ,'color', [color(subi,:)], 'LineWidth', 2);  %Residual
        plot(ppm, (zeros(1,length(ppm))), 'Color', 'k'); % Zeroline
        text(ppm(min_ppm_index)-0.05, 0, '0', 'FontSize', 18,'Color','k'); %Zeroline text

        %             plot(ppm,residual_mean+shift-maxshift_abs*0.3 ,'color', [color(subi,:)], 'LineWidth', 1);  %Residual
    end

%     if exist('baseline_mean', 'var')
%         plot(ppm,baseline_mean+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 2); %Baseline
%     end

    if exist('MM_mean', 'var')
        plot(ppm,MM_mean+baseline_mean+shift ,'color', colorFit, 'LineWidth', 2); %MM Baseline
    end

    if strcmp(which_spec,'mm') % re_mm
        if exist('MM_yu', 'var')
            fill([ppm fliplr(ppm)], [MM_yu+maxshift_abs*1.2 fliplr(MM_yl)+maxshift_abs*1.2], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow cleaned MM data
        end
        plot(ppm,MM_clean_mean+maxshift_abs*1.2 ,'color',[1 0 0.1], 'LineWidth', 3); % Data cleaned for MM data
    end



    %%% 6. DESIGN FINETUNING %%%
    % Adapt common style for all axes
    set(gca, 'XDir', 'reverse', 'XLim', [ppmRange(1), ppmRange(end)],'XMinorTick','on');
    ticks = get(gca,'XTick');
    set(gca, 'XTick', unique(round(ticks)));
    set(gca, 'LineWidth', 2, 'TickDir', 'out');
    set(gca, 'FontSize', 22);
    % If no y caption, remove y axis
    if isempty(ylab)
        set(gca, 'YColor', MRSCont.colormap.Background);
        set(gca,'YTickLabel',{});
        set(gca,'YTick',{});
        % Dirtywhite axes, light gray background
        set(gca, 'XColor', 'k');
        set(gca, 'Color', MRSCont.colormap.Background);
        set(gcf, 'Color', MRSCont.colormap.Background);
        %         title(figTitle, 'Interpreter', 'Tex', 'Color', 'k');
    else
        set(gca, 'YColor', 'k');
    end

    box off;
    xlabel(xlab, 'FontSize', 28);
    ylabel(ylab, 'FontSize', 28);
end
