function plot_moc(dir,sectionText,mocTop,mocLat,botDepth)

% Plot cross-sections of MPAS fields.
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% dir                text string, name of simulation
% sectionText        a cell array with text describing each section
% sectionID          section numbers for each row of this page
% pageName           name of this group of sections 
% mocData(nVertLevels,nPoints,nSections,nVars)
%   data in each cross-section for each variable
% depth(nVertLevels)                         depth of center of each layer, for plotting
% latSection(nPoints,nSections) lat coordinates of each section
% lonSection(nPoints,nSections) lon coordinates of each section
% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
% plotDepth(nSections) depth to which to plot each section
% var_lims(nVars,3)  contour line definition: min, max, interval 
% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% fid_latex           file ID of latex file

fprintf(['** sub_plot_cross_sections simulation: ' dir ...
	 '\n'])

px = [.54];
py=.53;
pw = [.87];  ph=[.83]; % width and height of plots
      temptext = char(sectionText(1));

% smooth once, like POP:
[nVertLevels nLat] = size(mocTop);
mocTopSmooth(:,1) = mocTop(:,1);
mocTopSmooth(:,nLat) = mocTop(:,nLat);
for j=2:nLat-1
  mocTopSmooth(:,j) = (mocTop(:,j-1) + mocTop(:,j-1) + mocTop(:,j-1))/3;
end

figure(3); clf

contour_lims = [-10:2:10];
        [cout,h]=contourf(mocLat,[0; botDepth],mocTopSmooth,...
               contour_lims);
	set(gca,'CLim',[min(contour_lims) max(contour_lims)])
	set(h,'LineColor',[.5 .5 .5])
        cbfrac=0;
	hold on

	% Text labels on countours
        [cout,h]=contour(mocLat,[0; botDepth],mocTopSmooth,...
			 contour_lims);
        ls=[200];
        clabel(cout,h,'fontsize',10,'color','k','rotation',0,'LabelSpacing',ls);
	set(h,'LineColor',[.5 .5 .5])

	% Black lines
        [cout,h]=contour(mocLat,[0; botDepth],mocTopSmooth,[-20:10:20]);
	set(h,'LineColor',[0 0 0],'LineWidth',1)

	
%contour(mocLat,[0 botDepth],mocTopSmooth,[-15:2:20])
set(gca,'YDir','reverse')
%colorbar
grid on
xlabel('latitude')
ylabel('depth')
%title(['Atlantic MOC, Sv, ' dir ', mean of years 11-20, solved from w'],'Interpreter','none');
title(['Atlantic MOC, Sv, ' dir ', 120km, solved from w'],'Interpreter','none');


	% stretched colorbar using contour_lims:
	cmin=min(contour_lims);
	cmax=max(contour_lims);
%	cvalue = cmin-.5*dc:dc:cmax+.5*dc;
        nc_orig = 256;
	nc = length(contour_lims);
	cmap_orig = ColdHot(nc_orig);
	cmap_orig_short = zeros(nc-1,3);
	ind=(.5:1:nc-.5);
	for j=1:nc-1
	  cmap_orig_short(j,:) = cmap_orig( floor((j-.5)/(nc-1)*nc_orig),:);
	end
	
	cvalue = linspace(cmin,cmax,256);
	nc_inc = length(cvalue);
	
	cmapnew = zeros(nc_inc,3);
	for jnew=2:nc_inc
	  jold = max(min(min(find(contour_lims>=cvalue(jnew))),nc)-1,1);
	  cmapnew(jnew-1,:) = cmap_orig_short(jold,:);
	end
	cmapnew(nc_inc,:) = cmap_orig_short(nc-1,:);
	
        colormap(cmapnew)

h=colorbar;
set(h,'YTick',[-16:4:16]);


   set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
     'PaperPosition',[0.25 0.25 7 3.2])
%   subplot('position',[0 .95 1 .05]); axis off
%   title_txt = [regexprep(char(var_name(iVar)),'_','\\_') ',  ' regexprep(dir,'_','\\_')];
%   h=text(.55,.4,title_txt);
%   set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%      text(.005,.7,[ date ]);

   unix(['mkdir -p f/' dir ]);
   temp=['f/' dir '/AtlanticMoc'];
   filename = regexprep(temp,'\.','_');
   print('-djpeg',[filename '.jpg']);
   print('-depsc2',[filename '.eps']);
   unix(['epstopdf ' filename '.eps --outfile=' filename '.pdf']);

%max(max(mocTopSmooth))
%min(min(mocTopSmooth))

return


px = [.28 .78];
py=linspace(.84,.13,4); % Midpoint position of plots
pw = [.4];  ph=[.18]; % width and height of plots

nPoints   = size(mocData,2);
nSections = size(mocData,3);

for iVar=1:1 %nVars
   %figure(iVar+1); clf
   %set(gcf,'Position',[100+(iVar*100) 600-iVar*100 550 400])
   %temptext2=char(var_name(iVar));

   for iRow = 1:length(sectionID)
      %figure(10 + iRow+1); clf
      %set(gcf,'Position',[100+(iRow*100) 1200-iRow*100 550 400])
      iSection = sectionID(iRow);
      if coord(iSection,1)==coord(iSection,3) % meridional section
	xtext = 'longitude';
	xaxis = lonSection(:,iSection);
      else % zonal section
	xtext = 'latitude';
	xaxis = latSection(:,iSection);
      end     

      % left column
  %    ha=subplot('position',[px(1)-pw/2 py(iRow)-ph/2 pw ph]);
  %    temptext = char(sectionText(iSection));
  %    if temptext2(1:6)=='ke_acc'
  %      h=surf(xaxis, depth,log10(mocData(:,1:nCellsInSection(iSection),iSection,iVar)));
  %      set(gca,'CLim',[-1 1.2])      
  %    else
  %      h=surf(xaxis, depth,mocData(:,1:nCellsInSection(iSection),iSection,iVar));
  %    end
     
  %    set(h,'EdgeColor','none')
  %    view(0,-90)
  %    title([temptext ', cm/s'])
  %    ylabel('depth, m')
  %    xlabel(xtext)
  %    axis tight
  %    set(gca,'YLim',[0 plotDepth(iSection)])
  %    h=colorbar  ;
  %    if temptext2(1:6)=='ke_acc'
%	set(h,'YTick',[-1:1:1.2],'YTickLabel',[0.1 1 10])
   %   end

      % right column
	
      hold on
%        contour(xaxis, depth,mocData(:,1:nCellsInSection(iSection),iSection,iVar), ...
%	       [var_lims(iVar,1):var_lims(iVar,3):var_lims(iVar,2)]);
%	set(gca,'CLim',var_lims(iVar,1:2))

        %%%%%% special commands for DWBC mrp
	% imitating colorbar at http://www.agu.org/journals/jc/jc1203/2011JC007586/figures.shtml#fig10
	
	%xaxis = xaxis - 360*ones(size(xaxis));
	% xaxis is now in longitude.  Convert to Distance (km)
        % along 26.5N east of 77W
	%xaxis = (xaxis+77)*99; % for DWBC only
	%contour_lims = [-25 -20 -15 -10 -5 -2 -1  1 2 5 10 15 20 25];
	%contour_lims = [-20 -15 -10 -5 -2 0 2 5 10 15 20 25 30];
	contour_lims = [var_lims(iVar,1):var_lims(iVar,3):var_lims(iVar,2)];
        [cout,h]=contourf(xaxis, depth,mocData(:,:,iSection,iVar), ...
               contour_lims);
	set(gca,'CLim',[min(contour_lims) max(contour_lims)])
	set(h,'LineColor',[.5 .5 .5])
        cbfrac=0;

	% Text labels on countours
        [cout,h]=contour(xaxis, depth,mocData(:,:,iSection,iVar),...
			 contour_lims);
        ls=[200];
        clabel(cout,h,'fontsize',10,'color','k','rotation',0,'LabelSpacing',ls);
	set(h,'LineColor',[.5 .5 .5])

	% Black lines
        [cout,h]=contour(xaxis, depth,mocData(:,:,iSection,iVar),[-100:50:100]);
	set(h,'LineColor',[0 0 0],'LineWidth',1)


	% stretched colorbar using contour_lims:
	cmin=min(contour_lims);
	cmax=max(contour_lims);
%	cvalue = cmin-.5*dc:dc:cmax+.5*dc;
        nc_orig = 256;
	nc = length(contour_lims);
	cmap_orig = ColdHot(nc_orig);
	cmap_orig_short = zeros(nc-1,3);
	ind=(.5:1:nc-.5);
	for j=1:nc-1
	  cmap_orig_short(j,:) = cmap_orig( floor((j-.5)/(nc-1)*nc_orig),:);
	end
	
	cvalue = linspace(cmin,cmax,256);
	nc_inc = length(cvalue);
	
	cmapnew = zeros(nc_inc,3);
	for jnew=2:nc_inc
	  jold = max(min(min(find(contour_lims>=cvalue(jnew))),nc)-1,1);
	  cmapnew(jnew-1,:) = cmap_orig_short(jold,:);
	end
	cmapnew(nc_inc,:) = cmap_orig_short(nc-1,:);
	
        colormap(cmapnew)
        %colorbarf_spec(cout,h,'vert',contour_lims);
        %xlabel('Distance (km) along 26.5N east of 77W') % for DWBC only

      axis tight
      %set(gca,'YLim',[0 plotDepth(iSection)],'XLim',[0 175]) % for DWBC only
      %set(gca,'YTick',[0:1000:5000],'XTick',[0:25:175])
      set(gca,'YLim',[0 plotDepth(iRow)])
      %set(gca,'YTick',[0:100:400])
      xlabel(xtext)
     % set(gca,'XTick',-1*[80:.5:70])
	%%%%%% special commands for DWBC mrp end

	%%%%%% special commands for EUC mrp end
        %if iRow==2
	%  set(gca,'XTick',[143 156 165 180 190   205   220   235   250   265])
	%  set(gca,'XTickLabel',{'143' '156' '165E' '180' '170W' '155' '140' '125' '110' '95'})
	%end
	
	%%%%%% special commands for EUC mrp end
	
        set(gca,'YDir','reverse')
        title([temptext ', ' num2str(var_lims(iVar,3)) ' Sv int'])
      ylabel('depth, m')
      grid on
      set(gca,'layer','top');
      h=colorbar;

      % mrp draw bottom based on zero contour
      hold on
      n = nPoints;
       % old way: maxLevelCell=zeros(1,n);
       x(2:n) = (xaxis(1:n-1)+xaxis(2:n))/2;
       x(1) = xaxis(1) - (xaxis(2)-xaxis(1))/2;
       x(n+1) = xaxis(n) + (xaxis(n)-xaxis(n-1))/2;
       b = max(botDepth);
       for j=1:n
          % old way: maxLevelCell(j)=max(min(find(mocData(:,j,iSection,iVar)==0.0))-1,1);
          %depthline(j) = botDepth(maxLevelCellSection(j,iSection));
	  %h=patch([x(j) x(j+1) x(j+1) x(j) x(j)],...
	  %        [b b depthline(j)  depthline(j) b], [.5 .5 .5]);
	  %set(h,'LineStyle','none')
       end
       
      % mrp draw bottom based on zero contour end

   set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
     'PaperPosition',[0.25 0.25 5.5 3.2])
   subplot('position',[0 .95 1 .05]); axis off
   title_txt = [regexprep(char(var_name(iVar)),'_','\\_') ',  ' regexprep(dir,'_','\\_')];
%   h=text(.55,.4,title_txt);
%   set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%      text(.005,.7,[ date ]);

   unix(['mkdir -p f/' dir ]);
   temp=['f/' dir '/' netcdfFile '_' pageName num2str(iRow) '_var' num2str(iVar)];
   filename = regexprep(temp,'\.','_');
   print('-djpeg',[filename '.jpg']);
   print('-depsc2',[filename '.eps']);
   unix(['epstopdf ' filename '.eps --outfile=' filename '.pdf']);
%   fprintf(fid_latex,['\\begin{figure}[btp]  \\center \n \\includegraphics[width=7.5in]{'...
%      filename '.jpg} \n\\end{figure} \n']);

   %   print('-depsc2',[filename '.eps'])
      
   end


end

