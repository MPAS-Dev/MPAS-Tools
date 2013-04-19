function sub_plot_section_locations(dir,coord, ...
       latSection,lonSection,fid_latex)

% Plot section locations on world map

% Mark Petersen, MPAS-Ocean Team, LANL, Sep 2012

%%%%%%%%%% input arguments %%%%%%%%%
% dir                text string, name of simulation
% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
% latSection(nPoints,nSections) lat coordinates of each section
% lonSection(nPoints,nSections) lon coordinates of each section
% fid_latex           file ID of latex file

fprintf(['** sub_plot_cell_sections, on figure 1.\n'])

nSections = size(coord,1);

figure(1); clf

  minLon = 0.0;
  latTrans = 360;
   
   % plot topo data of the earth.  This is just low-rez one deg
   % data for visual reference.
   load('topo.mat','topo','topomap1');
   if minLon==-180
     topoNew(:,1:180) = topo(:,181:360);
     topoNew(:,181:360) = topo(:,1:180);
     image([-180 180],[-90 90],topoNew,'CDataMapping', 'scaled');
   else
     image([0 360],[-90 90],topo,'CDataMapping', 'scaled');
   end     

   colormap(topomap1);
   set(gca,'YDir','normal')

   hold on

   % world
   axis tight
   set(gca,'XTick',30*[-10:12])
   set(gca,'YTick',15*[-20:20])

   % half world
   axis([-360+latTrans 0+latTrans -80 70]) 
   set(gca,'XTick',20*[-10:20])
   set(gca,'YTick',10*[-20:20])

   % N Atlantic
%   axis([-90+latTrans -5+latTrans -5 70]) 
%   set(gca,'XTick',[-100:5:360])
%   set(gca,'YTick',[-90:5:90])

   % Drake passage
%   axis([-90+latTrans,-50+latTrans,-75,-50]) 
%   set(gca,'XTick',[-100:2:360])
 %  set(gca,'YTick',[-200:2:200])

   % Pacific
%   axis([130 260 -10 10]) 
%   set(gca,'XTick',[0:1:300])
%   set(gca,'YTick',[-20:.1:20])
 
   hold on
   grid on

   for iSection=1:nSections
     h=plot(lonSection(:,iSection),latSection(:,iSection),'r-');
     h=text(lonSection(1,iSection),latSection(1,iSection), ...
	    num2str(iSection))
     get(h)
     
     set(h,'Color',[1 1 1],'FontWeight','bold')
     %h=plot(lonSection(:,iSection),latSection(:,iSection),'y.');
     %set(h,'Color','y','LineWidth',1)
   end

   ylabel('latitude')
   xlabel('longitude')
   title(['Domain: ' regexprep(dir,'_','\\_') ' Cells in cross sections. '])

   set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	   'PaperPosition',[0.25 0.25 8 8])

   subplot('position',[0 .95 1 .05]); axis off
   text(.005,.7,[ date ]);

   dir_name1 =  regexprep(dir,'\.','_');
   dir_name2 =  regexprep(dir_name1,'/','_');
   filename=['f/' dir_name2 '_cell_map' ];
   print('-djpeg',[filename '.jpg'])
   
   % put printing text in a latex file
   fprintf(fid_latex,...
     ['\\begin{figure}[btp]  \\center \n \\includegraphics[width=7.5in]{'...
      filename '.jpg} \n\\end{figure} \n']);
