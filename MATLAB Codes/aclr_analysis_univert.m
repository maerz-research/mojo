% ACLR_ANALYSIS_UNIVERT.M is a script which analyzes mouse ACL rupture 
% mechanical data collected on the CellScale Univert system. The correct
% folder structure for this script is as follows:
% 
% - Root Folder
%     - ACLR Mechanical Data
%         - Subfolders containing CellScale outputs for each test - these
%           can be dragged and dropped from the CellScale Data folder
%           after each session of ruptures
%     - MATLAB Analysis (set as working directory)
%         - (This script)
%         - Graphs
%             - This folder is where visual outputs from this script will
%               print to
% 
% If folder/file organization is correct, this script should automatically
% analyze all available datasets (skipping any datasets which have been
% previously analyzed) and output the following:
% 
% 1. Graphs showing the ACL rupture event, preload/preconditioning, and ACL
%    rupture steady-state stiffness (i.e. the linear region of the force/
%    displacement curve)
% 2. Numerical parameters about the test, including failure load,
%    stiffness, etc. These are aggregated into running .MAT and .CSV files
% 
% Version history:
% v1.0 - Jan 2024 - Basic functionality for analyzing ACL rupture event in
%                   a loop.
% v1.1 - Apr 2024 - Added a check to exclude already-analyzed datapoints, 
%                   additional analysis of preload/preconditioning phases 
%                   with an additional graph, and changed data storage from
%                   cell to table.
% 
% Contact: Mike Newton, University of Michigan, nemichae@umich.edu.

%% Initialize paths and directories

% Create directory of data - data should be stored
p1 = [pwd,'\..\ACLR Mechanical Data\'] ; % Data directory path
d  =  dir2cell(dir([p1,'*ACLR*']),true); % Directory of samples

% Ensure graph directories are already present
if ~exist([pwd,'\Graphs'],'dir'), mkdir([pwd,'\Graphs']); end

% Load in running results compilation, formatted as a table
if exist('compiled_results.mat','file')
    r = load('compiled_results.mat'); results_all = r.results_all; clear r;
else, results_all = {'Name','Failure Load (N)','Failure Displacement (mm)',...
                     'Total Displacement (mm)','Stiffness (N/mm)',...
                     'Steady-State Velocity (mm/s)','Total Creep (mm)',...
                     '% Error - Precond 1','% Error - Cyc Min',...
                     '% Error - Cyc Max','% Error - Precond 2'};
end
if iscell(results_all)
    results_all = cell2table(results_all(2:end,:),'VariableNames',results_all(1,:));
end

d = d(~ismember(d,results_all.Name)); % Remove samples already processed
L = length(d);

%%
for ii = 1:L
    %%
    t = tic;

    % Read in data for current run
    d2 = dir2cell(dir([p1,d{ii},'\*.csv']));
    if isempty(d2), continue; end
    tab = readtable([p1,d{ii},'\',d2{1}],'VariableNamingRule','preserve');
    tab = tab(:,1:6);

    % Section 1: Outcomes of ACL rupture event

    % Find the start and end of the rupture displacement
    I_aclr =  find(ismember(tab.SetName,'Downward-Displacement') & ...
                   ismember(tab.Cycle  ,'1-Compress'           ));
    r_aclr = [min(I_aclr),max(I_aclr)];
    if isempty(r_aclr), continue; end

    % Account for the 1-Recover glitch, where the beginning of the
    % displacement begins within the end of the preload in an added
    % '1-Recover' section rather than on the first line of the
    % 'Downward-Displacement' section
    if strcmp(tab.Cycle{r_aclr(1)-1},'1-Recover')
        r_aclr(1) = find(ismember(tab.SetName,'Preload-1N-2') & ...
                         ismember(tab.Cycle  ,'1-Recover'   ),1);
    end

    % Find the starting and end of the retraction for graphing
    I_retr =  find(ismember(tab.SetName,'Retract') & ...
                   tab.Displacement_mm >= tab.Displacement_mm(r_aclr(1)));
    r_retr = [min(I_retr),max(I_retr)];

    time  = tab.Time_S(         r_aclr(1):r_retr(2))-tab.Time_S(         r_aclr(1));
    disp  = tab.Displacement_mm(r_aclr(1):r_retr(2))-tab.Displacement_mm(r_aclr(1));
    force = tab.Force_N(        r_aclr(1):r_retr(2))                               ;

%     % Find the first zero-crossing of the velocity curve - this should
%     % correspond to the frame on which ACL rupture occurs
%     veloc = diff(medfilt1(force),5)      ;
%     fail  = find(diff(sign(veloc)),1) + 2; % +2 accounts for 2 diff() fxns
%     fail  = fail - 3 + find(force(fail-2:fail+2)==max(force(fail-2:fail+2))); % Fine-tune

    % Use MATLAB's 'findchangepts' function (Signal Processing Toolbox) to
    % identify where the rupture occurred - the function *should* assign
    % the failure to the rupture when one occurs, and to the final change
    % in displacement when a rupture does not occur. Because occasionally
    % there is a ramp-up acceleration phase in the first few frames of the
    % displacement, the first 6 frames are skipped because they can
    % occasionally be false-positively identified as a "change-point" and
    % thereby affect the analysis
    fail = findchangepts(force(1:diff(r_aclr)+5),'Statistic','Linear','MaxNumChanges',1);
    fail = fail(1) - 6 + find(force(fail(1)-5:fail(1)+5)==...
                          max(force(fail(1)-5:fail(1)+5))); % Fine-tune
    fail = fail(1);

    % Perform stiffness analysis middle 50% of start-to-failure data
    bumper      = round(fail * 0.25);
    disp_stiff  = disp( 1+bumper:fail-bumper);
    force_stiff = force(1+bumper:fail-bumper);
    time_stiff  = time( 1+bumper:fail-bumper);
    linfit      = polyfit( disp_stiff,force_stiff,1);
    R2          = corrcoef(disp_stiff,force_stiff  ); R2 = R2(2)^2;
    linfit2     = polyfit( time_stiff, disp_stiff,1);

    % Visualize stiffness calculation
    figure,scatter(disp_stiff,force_stiff,'filled','MarkerFaceColor',[0,.6,0]); hold on
    plot(disp_stiff,polyval(linfit,disp_stiff),'Color','k','LineWidth',1,...
        'LineStyle','--');
    xlabel('Displacement (mm)'); ylabel('Load (N)');
    dim = [0.2 0.5 0.3 0.3];
    str = {['y = ' ,num2str(linfit(1),'%.2f'),'x + ',num2str(linfit(2),'%.2f')],...
           ['R2 = ',num2str(R2,'%.3f')]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
    print(['Graphs\stiffness_',d{ii},'.tif'],'-dtiff','-r150');
    close(gcf);

    % Visualize test section
    figure;
    subplot(2,1,1); plot(time,disp ,'LineWidth',1,'Color','b');
    xlabel('Time (s)'); ylabel('Displacement (mm)');
    ylim([0,2]); hold on; plot(get(gca,'xlim'),[1.5,1.5],...
        'LineWidth',1,'Color','k','LineStyle','--');
    scatter(time(fail),disp(fail),'*','MarkerEdgeColor','b');
    subplot(2,1,2); plot(time,force,'LineWidth',1,'Color','r');
    xlabel('Time (s)'); ylabel('Load (N)');
    hold on; plot([0.15,0.15],get(gca,'ylim'),...
        'LineWidth',1,'Color','k','LineStyle','--');
    scatter(time(fail),force(fail),'*','MarkerEdgeColor','r');
    print(['Graphs\aclr_',d{ii},'.tif'],'-dtiff','-r150');
    close(gcf);

    % Derive outcome measures
    fail_load  = force(fail)                      ;
    fail_disp  = disp( fail)                      ;
    total_disp = max(tab.Displacement_mm(r_aclr(1):r_aclr(2))) - ...
                     tab.Displacement_mm(r_aclr(1));
    stiffness  = linfit( 1)                       ;
    ss_veloc   = linfit2(1)                       ;

    % Section 2: Characterize Preload and Preconditioning

    % Find the start and end of the rupture displacement
    I_pre1   =  find(ismember(tab.SetName,'Preload-1N'  ) & ...
                     ismember(tab.Cycle  ,'1-Hold'));
    r_pre1   = [min(I_pre1),max(I_pre1)];
    I_pre2   =  find(ismember(tab.SetName,'Preload-1N-2') & ...
                     ismember(tab.Cycle  ,'1-Hold'));
    r_pre2   = [min(I_pre2),max(I_pre2)];
    r_cyc    = zeros(10,4);
    cyc_mnmx = zeros(10,4);
    for jj = 1:10
        I_cyc1 = find(ismember(tab.SetName,'1-3N Cycle') & ...
                      ismember(tab.Cycle  ,[num2str(jj),'-Compress']));
        I_cyc2 = find(ismember(tab.SetName,'1-3N Cycle') & ...
                      ismember(tab.Cycle  ,[num2str(jj),'-Recover']));
        r_cyc(   jj,:) = [min(I_cyc1),max(I_cyc1),min(I_cyc2),max(I_cyc2)];
        cyc_mnmx(jj,:) = [min(tab.Force_N(        r_cyc(jj,3):r_cyc(jj,4))),...
                          max(tab.Force_N(        r_cyc(jj,1):r_cyc(jj,2))),...
                          min(tab.Displacement_mm(r_cyc(jj,3):r_cyc(jj,4))),...
                          max(tab.Displacement_mm(r_cyc(jj,1):r_cyc(jj,2)))];
    end
    cyc_mnmx(:,5) = mean(cyc_mnmx(:,3:4),2);

    % Visualize test section
    figure;
    subplot(2,1,1); plot(tab.Time_S(         r_pre1(1):r_pre2(end)),...
                         tab.Displacement_mm(r_pre1(1):r_pre2(end)),...
                        'LineWidth',1,'Color','b'); hold on;
    plot(tab.Time_S(round(mean([r_cyc(:,1),r_cyc(:,4)],2))),...
         cyc_mnmx(:,5),'LineWidth',1,'LineStyle','--','Color','k');
    plot(tab.Time_S(         [r_pre1(1),r_pre1(2)]),...
         tab.Displacement_mm([r_pre1(1),r_pre1(1)]),...
        'LineWidth',1,'LineStyle','--','Color','k');
    plot(tab.Time_S(         [r_pre2(1),r_pre2(2)]),...
         tab.Displacement_mm([r_pre2(2),r_pre2(2)]),...
        'LineWidth',1,'LineStyle','--','Color','k');
    xlabel('Time (sec)'); ylabel('Displacement (mm)');
    subplot(2,1,2); plot(tab.Time_S( r_pre1(1):r_pre2(end)),...
                         tab.Force_N(r_pre1(1):r_pre2(end)),...
                        'LineWidth',1,'Color','r'); hold on;
    plot(tab.Time_S([r_cyc(1),r_cyc(end)]),[3,3],...
        'LineWidth',1,'LineStyle','--','Color','k');
    plot(tab.Time_S([r_pre1(1),r_pre2(end)]),[1,1],...
        'LineWidth',1,'LineStyle','--','Color','k');
    ylim([0,4]); xlabel('Time (sec)'); ylabel('Load (N)');
    print(['Graphs\precond_',d{ii},'.tif'],'-dtiff','-r150');
    close(gcf);

    % Derive outcome measures
    total_creep = tab.Displacement_mm(r_pre2(end))-tab.Displacement_mm(r_pre1(1));
    % total_creep = mean(tab.Displacement_mm(r_pre2(2)-100:r_pre2(2)))...
    %                        - tab.Displacement_mm(r_pre1(1));
    pe_pre1 = mean(abs(1-(tab.Force_N(r_pre1(1):r_pre1(2)))));
    pe_pre2 = mean(abs(1-(tab.Force_N(r_pre2(1):r_pre2(2)))));
    pe_min  = mean(abs(1-cyc_mnmx(:,1)))  ;
    pe_max  = mean(abs(3-cyc_mnmx(:,2)))/3;

    % Add current run's results to compiled results
    results_current = cell2table({d{ii},fail_load,fail_disp,total_disp,...
                             stiffness,ss_veloc,total_creep,pe_pre1,...
                             pe_min,pe_max,pe_pre2},'VariableNames',...
                             results_all.Properties.VariableNames);
    results_all = [results_all;results_current];

    looptrack(ii,L,t,d{ii});
end

% Save results cell
results_all = sortrows(unique(results_all,'rows'),1);
last_update = date;
save('compiled_results.mat','results_all','last_update');
writetable(results_all,'compiled_results.csv');