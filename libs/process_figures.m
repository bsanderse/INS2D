%% some more helpful commands:
% colors of the standard Matlab colors:
% initialColorOrder = get(gca,'ColorOrder'); 
% h = figure(1) % or whatever fig
% get the axes handle
% a = get(h, 'CurrentAxes')
% get the handles for children of the axes -- these are the data series handles
% c = get(a, 'Children')         
% set(c(1),'Color',initialColorOrder(4,:))

% set font size:
% set(gca,'FontSize',13);

% set background color to white
% set(gcf,'color','w');


%%
% load all existing figures in folder and increase marker size (or other property)

close all
clear all

folder = uigetdir; % check the help for uigetdir to see how to specify a starting path, which makes your life easier

% get the names of all files. dirListing is a struct array. 
dirListing = dir(folder);

% loop through the files and open. Note that dir also lists the directories, so you have to check for them.
for d = 1:length(dirListing)
    
    if (~dirListing(d).isdir)
        fileName = fullfile(folder,dirListing(d).name); % use full path because the folder may not be the active path
        
        if (strcmp(fileName(end-2:end),'fig') & strfind(fileName,'nocorr'))
            % open your file here 
            open(fileName)
            %%
            h = figure(1);
            a = get(h, 'CurrentAxes');
            % get the handles for children of the axes -- these are the data series handles
            c = get(a, 'Children');            
            % this is strange, but we have to do this twice to make it work
            a = get(h, 'CurrentAxes');
            c = get(a, 'Children');
            
            cmap = get(gca,'ColorOrder');

            n = length(c);
            styles= {'x','--','-.',':','-','--'};
%             color = {'k','k','g','c','r','b'};

            for i=1:n

%                 set(c(n+1-i),'MarkerSize',10);
                set(c(n+1-i),'LineWidth',1);
%                 set(c(n+1-i),'LineStyle',styles{i});
%                 set(c(n+1-i),'Color',color{i});
                set(c(n+1-i),'Color',cmap(i,:)); 
                
            end
            set(c(n),'Marker','s')
            layout(1,18,'p','x');
            
%             set(gca,'YTick',[1e-14 1e-10 1e-6 1e-2]);

%             if (strfind(fileName,'energy_error'))
%                 disp(fileName)
%                 set(gca,'YTick',[1e-16 1e-12 1e-8 1e-4]);
%             end
%             if (strfind(fileName,'timerev_error'))
%                 disp(fileName)
%                 set(gca,'YTick',[1e-14 1e-10 1e-6 1e-2]);
%             end
%             if (strfind(fileName,'enstrophy_error_t8'))
%                 disp(fileName)
%                 set(gca,'YTick',[1e-5 1e-3 1e-1]);                
%             end
            %%
            file_eps = [fileName(1:end-3) 'eps'];
            % save as eps
            print(gcf,'-depsc',file_eps);
            % save as fig
%             saveas(gcf,fileName,'fig');
            close(1)
            
        end
    end
    
end 
