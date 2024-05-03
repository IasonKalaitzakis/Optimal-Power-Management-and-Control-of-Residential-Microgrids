function SubPlot
% -----------------------------------------------------------
% General usage: This function inserts various MATLAB figure (.fig) files
%                into one figure with multiple subplots
%
% Notes:
%      All desired .fig files should contain only one figure inside it (no
%      subplot inside the .fig files).
%      All desired .fig files should be defined in the 2D space.
%
% Author: Farhad Sedaghati
% Contact: <farhad_seda@yahoo.com>
% Written: 04/08/2015
% Updated: 06/21/2015
% Updated: 07/14/2015
% Updated: 08/13/2015
% Updated: 10/08/2019 to include the multiselect option, and get the legend,
% x and y limit, scale, tick, and ticklabels
% -----------------------------------------------------------


disp('---------------------------------------------------------------------');
disp('note that all desired .fig files should contain only one figure inside it.');
disp('note that all desired .fig files should be defined in the 2D space.');
disp('Please press enter to continue');
pause


% Number of the run
N=0;
answer='y';
while strcmpi(answer,'y')
    
    % Get the path and filename of the desired fig file
    filename=0;
    
    run=0;
    while isequal(filename,0)
        if run==2
            error('please choose your fig file');
        end
        disp('---------------------------------------------------');
        disp('Select the desired fig file(files) which you want to insert in the subplot: ');
        [filename,pathname]=uigetfile('MultiSelect', 'on',{'*.fig';'*.FIG'},'Select the .fig file(files) you want to insert in the subplot');
        run=run+1;
    end
    clear FN
    
    if ~iscell(filename)
        N=N+1;
        FN=[pathname filename];
        % open figure
        h(N) = openfig(FN,'new');
        % get handle to axes of figure
        ax(N)=gca;
    else
        for figureIndex = 1:size(filename,2)
            N=N+1;
            FN=[pathname filename{1,figureIndex}];
            % open figure
            h(N) = openfig(FN,'new');
            % get handle to axes of figure
            ax(N)=gca;
        end
    end
    answer=input('Do you have more .fig files to read? \n','s');
end


K=input('How many figures in a row do you want to have? \n');
if isempty(K)
    K=2;
end


figure;
for i=1:N
    % create and get handle to the subplot axes
    s(i) = subplot(ceil(N/K),K,i); 
    % get handle to all the children in the figure
    aux=get(ax(i),'children');
    for j=1:size(aux)
        fig(i) = aux(j);
        copyobj(fig(i),s(i)); 
        hold on
    end
    % copy children to new parent axes i.e. the subplot axes
    xlab = get(get(ax(i),'xlabel'),'string');
    ylab = get(get(ax(i),'ylabel'),'string');
    tit = get(get(ax(i),'title'),'string');
    Legend = get(get(ax(i),'legend'),'string');
    xLimits = get(ax(i),'XLim');
    yLimits = get(ax(i),'YLim');
    ScaleX = get(ax(i),'Xscale');
    ScaleY = get(ax(i),'Yscale');
    TickX = get(ax(i),'Xtick');
    TickY = get(ax(i),'Ytick');
    TicklabelX = get(ax(i),'Xticklabel');
    TicklabelY = get(ax(i),'Yticklabel');
    set(gca, 'XScale', ScaleX, 'YScale', ScaleY, 'Xtick', TickX, 'Ytick',...
        TickY, 'Xticklabel', TicklabelX, 'Yticklabel', TicklabelY);
    xlabel(xlab);ylabel(ylab);title(tit);xlim(xLimits);ylim(yLimits);
    legend(Legend)
end








