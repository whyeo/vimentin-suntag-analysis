function formatfig(fig,style,savename,painters)

savefile = true;

if nargin <= 0
    fig = gcf;
end

if nargin <= 1
    style = 'presentation';
end

if nargin <= 2
    savefile = false;
elseif ~isstring(savename)
    if islogical(savename)
        savefile = savename;
    end
end

if nargin <= 3
    painters = false;
elseif ischar(painters) && strcmpi(painters,'painters')
    painters = true;
end

style_list = {'presentation', 'publication'};
assert(contains(style,style_list));

ax = findobj(fig, 'type', 'axes');

try
    for ii = 1:length(ax)
        if contains(style, 'presentation')
            if ismac
                ax(ii).FontSize = 16;
            else
                ax(ii).FontSize = 12;
            end
        elseif contains(style, 'publication')
            if ismac
                ax(ii).FontSize = 10;
            else
                ax(ii).FontSize = 8;
            end
        end
    end
    
catch
    disp('reverting to default method for fontsize');
    if contains(style, 'presentation')
        if ismac
            ax(ii).FontSize = 16;
        else
            ax(ii).FontSize = 12;
        end
    elseif contains(style, 'publication')
        if ismac
            ax(ii).FontSize = 10;
        else
            ax(ii).FontSize = 8;
        end
    end
end

if contains(style, 'presentation')
    set(fig,'DefaultLineLineWidth',1.5);
elseif contains(style, 'publication')
    set(fig,'DefaultLineLineWidth',1);
end

if savefile
    if endsWith(savename,'png')
        saveas(fig,savename);
    elseif endsWith(savename,'eps') || endsWith(savename,'epsc')
        if painters
            set(fig,'Renderer','Painters');
        end
        saveas(fig,savename,'epsc');
    elseif endsWith(savename,'svg')
        if painters
            set(fig,'Renderer','Painters');
        end
        saveas(fig,savename,'svg');
    end
end