function sub_plot_edge_sections(dir,sectionCoord, ...
     latSectionVertex,lonSectionVertex, ...
     latVertexDeg,lonVertexDeg, ...
     sectionEdgeIndex, nEdgesInSection,...
     fid_latex)

% Plot edge section locations on world map

% Mark Petersen, MPAS-Ocean Team, LANL, May 2012

%%%%%%%%%% input arguments %%%%%%%%%
% dir                text string, name of simulation
% sectionCoord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
% latVertexDeg(nVertices)                         lat arrays for all vertices
% lonVertexDeg(nVertices)                         lon arrays for all vertices  
% sectionEdgeIndex(maxEdges,nSections)       edge index of each section
% nEdgesInSection(nSections)                 number of edges in each section
% fid_latex           file ID of latex file

fprintf(['** sub_plot_edge_sections, on figure 1.\n'])

nSections = size(sectionCoord,1);

figure(1); clf

if (min(lonVertexDeg)>-1e-8)
  minLon = 0.0;
  latTrans = 360;
else
  minLon = -180.0;
  latTrans = 0.0;
end
   
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
   patch([-10 1000 1000 -10 -10],[-100 -100 100 100 -100],[.5 1 0])
   patch([-10 1000 1000 -10 -10],[-100 -100 100 100 -100],[1 1 1])
   set(gca,'YDir','normal')

   hold on

   % world
   axis([0 360 -90 90])
   set(gca,'XTick',30*[-10:12])
   set(gca,'YTick',15*[-20:20])

   % half world
%   axis([-240+latTrans 0+latTrans -80 70]) 
%   set(gca,'XTick',20*[-10:20])
%   set(gca,'YTick',10*[-20:20])

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
 
 
   % plot vertexs.  This is just done for debugging.
   h=plot(lonVertexDeg,latVertexDeg,'.b');
   set(h,'MarkerSize',2)

   grid on

   for iSection=1:nSections
     latCoordDeg = [sectionCoord(iSection,1) sectionCoord(iSection,3)];
     lonCoordDeg = [sectionCoord(iSection,2) sectionCoord(iSection,4)];

     %h=plot([mod(lonCoordDeg,360)],[latCoordDeg],'*-');
     %set(h,'Color','y','LineWidth',1)
     %h=plot([mod(lonCoordDeg(1),360)],[latCoordDeg(1)],'*k');

     for i=1:nEdgesInSection(iSection)
	h = line([lonSectionVertex(i,iSection) lonSectionVertex(i+1,iSection)],...
		 [latSectionVertex(i,iSection) latSectionVertex(i+1,iSection)]);
	set(h,'Color','r','LineWidth',2)
	%plot([lonVertexDeg(sectionVertexIndex(i+1,iSection))], ...
	%     [latVertexDeg(sectionVertexIndex(i+1,iSection))],'sk')
     end
   end

   ylabel('latitude')
   xlabel('longitude')
   title(['Domain: ' regexprep(dir,'_','\\_') ' Edges of transport sections. '])

   set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	   'PaperPosition',[0.25 0.25 16 8])

   subplot('position',[0 .95 1 .05]); axis off
   text(.005,.7,[ date ]);

   dir_name1 =  regexprep(dir,'\.','_');
   dir_name2 =  regexprep(dir_name1,'/','_');
   filename=['f/' dir_name2 '_vertex_map' ];
   print('-djpeg',[filename '.jpg']);
   
   % put printing text in a latex file
   fprintf(fid_latex,...
     ['\\begin{figure}[btp]  \\center \n \\includegraphics[width=7.5in]{'...
      filename '.jpg} \n\\end{figure} \n']);
