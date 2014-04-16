function plot_moc(dir,sectionText,mocTop,mocLat,botDepth, ...
		  contour_lims,var_name, fign)

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
% contour_lims(nVars,3)  contour line definition: min, max, interval 
% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% fid_latex           file ID of latex file

fprintf(['** sub_plot_cross_sections simulation: ' dir ...
	 '\n'])

px = [.54];
py=.53;
pw = [.87];  ph=[.83]; % width and height of plots

% smooth once, like POP:
[nVertLevels nLat] = size(mocTop);
mocTopSmooth(:,1) = mocTop(:,1);
mocTopSmooth(:,nLat) = mocTop(:,nLat);
for j=2:nLat-1
  mocTopSmooth(:,j) = (mocTop(:,j-1) + mocTop(:,j-1) + mocTop(:,j-1))/3;
end

figure(fign); clf

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
        %[cout,h]=contour(mocLat,[0; botDepth],mocTopSmooth,[-100:100:100]);
	%set(h,'LineColor',[0 0 0],'LineWidth',1)

	
%contour(mocLat,[0 botDepth],mocTopSmooth,[-15:2:20])
set(gca,'YDir','reverse')
%colorbar
grid on
xlabel('latitude')
ylabel('depth')
title([char(sectionText(1)) ', Sv, ' dir ', ' char(var_name(1))],'Interpreter','none');


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
set(h,'YTick',contour_lims);


   set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
     'PaperPosition',[0.25 0.25 7 3.2])
%   subplot('position',[0 .95 1 .05]); axis off
%   title_txt = [regexprep(char(var_name(iVar)),'_','\\_') ',  ' regexprep(dir,'_','\\_')];
%   h=text(.55,.4,title_txt);
%   set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%      text(.005,.7,[ date ]);

   unix(['mkdir -p f/' dir ]);
   tempTxt = char(sectionText(1));
   temp=['f/' dir '/' tempTxt(1:6) 'Moc_' char(var_name(1))];
   filename = regexprep(temp,'\.','_');
   print('-djpeg',[filename '.jpg']);
   print('-depsc2',[filename '.eps']);
   unix(['epstopdf ' filename '.eps --outfile=' filename '.pdf']);
