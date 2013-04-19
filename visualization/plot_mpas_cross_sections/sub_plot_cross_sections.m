function sub_plot_cross_sections(dir,sectionText,pageName,sectionID,sectionData,depth,...
   sectionCellIndex, nCellsInSection, latSection,lonSection, coord, plotDepth,...
   var_name,var_lims,fid_latex)

% Plot cross-sections of MPAS fields.
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% dir                text string, name of simulation
% sectionText        a cell array with text describing each section
% sectionID          section numbers for each row of this page
% pageName           name of this group of sections 
% sectionData(nVertLevels,max(nCellsInSection),nSections,nVars)
%   data in each cross-section for each variable
% depth(nVertLevels)                         depth of center of each layer, for plotting
% sectionCellIndex(maxCells,nSections)       cell index of each section
% nCellsInSection(nSections)                 number of cells in each section
% latSection(max(nCellsInSection),nSections) lat coordinates of each section
% lonSection(max(nCellsInSection),nSections) lon coordinates of each section
% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
% plotDepth(nSections) depth to which to plot each section
% var_lims(nVars,3)  contour line definition: min, max, interval 
% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% fid_latex           file ID of latex file

fprintf(['** sub_plot_cross_sections simulation: ' dir ...
	 ' plotting page: ' pageName '\n'])

px = [.28 .78];
py=linspace(.84,.13,4); % Midpoint position of plots
pw = [.4];  ph=[.18]; % width and height of plots
nSections = length(nCellsInSection);
nVars = length(var_name);

for iVar=1:nVars
   figure(iVar+1); clf
   set(gcf,'Position',[100+(iVar*100) 600-iVar*100   715   975])
   temptext2=char(var_name(iVar));

   for iRow = 1:length(sectionID)
      iSection = sectionID(iRow);
      if coord(iSection,1)==coord(iSection,3) % meridional section
	xtext = 'longitude';
	xaxis = lonSection(1:nCellsInSection(iSection),iSection);
      else % zonal section
	xtext = 'latitude';
	xaxis = latSection(1:nCellsInSection(iSection),iSection);
      end     

      % left column
      ha=subplot('position',[px(1)-pw/2 py(iRow)-ph/2 pw ph]);
      temptext = char(sectionText(iSection));
      if temptext2(1:6)=='ke_acc'
        h=surf(xaxis, depth,log10(sectionData(:,1:nCellsInSection(iSection),iSection,iVar)));
        set(gca,'CLim',[-1 1.2])      
      else
        h=surf(xaxis, depth,sectionData(:,1:nCellsInSection(iSection),iSection,iVar));
      end
     
      set(h,'EdgeColor','none')
      view(0,-90)
      title([temptext ', cm/s'])
      ylabel('depth, m')
      xlabel(xtext)
      axis tight
      set(gca,'YLim',[0 plotDepth(iSection)])
      h=colorbar  ;
      if temptext2(1:6)=='ke_acc'
	set(h,'YTick',[-1:1:1.2],'YTickLabel',[0.1 1 10])
      end

      % right column
	
      ha=subplot('position',[px(2)-pw/2 py(iRow)-ph/2 pw ph]);
      temptext = char(sectionText(iSection));
      hold on
        contour(xaxis, depth,sectionData(:,1:nCellsInSection(iSection),iSection,iVar), ...
	       [-30:var_lims(iVar,3):30]);
        set(gca,'CLim',var_lims(iVar,1:2))
        set(gca,'YDir','reverse')
        title([temptext ', ' num2str(var_lims(iVar,3)) 'cm/s int'])
      ylabel('depth, m')
      xlabel(xtext)
      axis tight
      set(gca,'YLim',[0 plotDepth(iSection)])
      grid on
      h=colorbar;
   
   end

   set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
     'PaperPosition',[0.25 0.25 8 10])
   subplot('position',[0 .95 1 .05]); axis off
   title_txt = [regexprep(char(var_name(iVar)),'_','\\_') ',  ' regexprep(dir,'_','\\_')];
   h=text(.55,.4,title_txt);
   set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
      text(.005,.7,[ date ]);

   dir_name1 =  regexprep(dir,'\.','_');
   dir_name2 =  regexprep(dir_name1,'/','_');
   filename=['f/' dir_name2 '_' pageName '_var' num2str(iVar)];
   print('-djpeg',[filename '.jpg'])
   fprintf(fid_latex,['\\begin{figure}[btp]  \\center \n \\includegraphics[width=7.5in]{'...
      filename '.jpg} \n\\end{figure} \n']);

   %   print('-depsc2',[filename '.eps'])

end

