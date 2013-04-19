function write_edge_sections_text ...
     (dir, ...
      sectionEdgeIndex, sectionEdgeSign, nEdgesInSection,...
      sectionText,sectionAbbreviation,sectionCoord)

% Write section edge information to the text file
% text_files/your_dir_transport_section_edges.nc
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% dir                string with run directory name
% sectionEdgeIndex(maxNEdgesInSection,nSections) edge index of each section
% sectionEdgeSign(maxNEdgesInSection,nSections)  sign of each
%    section, positive is to right of path direction.
% nEdgesInSection(nSections)                 number of cells in each section
% sectionText        a cell array with text describing each section
% sectionAbbreviation an 8-character title for each section
% sectionCoord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Write section edge information to a text file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
fprintf(['** Write edge information to file: ' dir '\n'])

nSections = length(nEdgesInSection);
maxNEdgesInSection = max(nEdgesInSection);

dir_name1 =  regexprep(dir,'/','_');
unix(['mkdir -p text_files/' dir_name1 ]);

% sectionEdgeIndex
filename = ['text_files/' dir_name1 '/sectionEdgeIndex.txt'];
fid = fopen(filename,'w');
for j=1:nSections
  fprintf(fid,' %10i',sectionEdgeIndex(:,j));
  fprintf(fid,' \n');
end
fclose(fid);

% sectionEdgeIndex
dir_name1 =  regexprep(dir,'/','_');
filename = ['text_files/' dir_name1 '/sectionEdgeIndex.txt'];
fid = fopen(filename,'w');
for j=1:nSections
  fprintf(fid,' %10i',sectionEdgeIndex(:,j));
  fprintf(fid,' \n');
end
fclose(fid);

% sectionEdgeSign
dir_name1 =  regexprep(dir,'/','_');
filename = ['text_files/' dir_name1 '/sectionEdgeSign.txt'];
fid = fopen(filename,'w');
for j=1:nSections
  fprintf(fid,' %10i',sectionEdgeSign(:,j));
  fprintf(fid,' \n');
end
fclose(fid);

% nEdgesInSection
dir_name1 =  regexprep(dir,'/','_');
filename = ['text_files/' dir_name1 '/nEdgesInSection.txt'];
fid = fopen(filename,'w');
for j=1:nSections
  fprintf(fid,' %10i',nEdgesInSection(j));
  fprintf(fid,' \n');
end
fclose(fid);

% sectionText
dir_name1 =  regexprep(dir,'/','_');
filename = ['text_files/' dir_name1 '/sectionText.txt'];
fid = fopen(filename,'w');
for j=1:nSections
  fprintf(fid,' %s',char(sectionText(j)));
  fprintf(fid,' \n');
end
fclose(fid);

% sectionAbbreviation
dir_name1 =  regexprep(dir,'/','_');
filename = ['text_files/' dir_name1 '/sectionAbbreviation.txt'];
fid = fopen(filename,'w');
for j=1:nSections
  fprintf(fid,' %s',sectionAbbreviation(j,:));
  fprintf(fid,' \n');
end
fclose(fid);

% sectionCoord
dir_name1 =  regexprep(dir,'/','_');
filename = ['text_files/' dir_name1 '/sectionCoord.txt'];
fid = fopen(filename,'w');
for j=1:nSections
  fprintf(fid,' %10.3f',sectionCoord(j,:));
  fprintf(fid,' \n');
end
fclose(fid);

