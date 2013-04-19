% read_u_TransportSection
% 
% Matlab script to plot output of uTransportSection.txt,
% the high frequency instantanious transport output.

% Mark Petersen, MPAS-Ocean Team, LANL, May 2012

wd = '.';
filename=[wd '/uTransportSection.txt'];

a= load('uTransportSection.txt');
nTimeSlices = size(a,1);
nSections = size(a,2)-1;
uTransportSection = a(:,2:nSections+1);

% this is not quite working.

%sectionAbbreviation = load(filename);
% sectionAbbreviation
%filename=[wd '/sectionAbbreviation.txt'];
%fid = fopen(filename,'r');
%for j=1:nSections
%  for k=1:8
%     sectionAbbreviation(j,k) = fscanf(fid,' %c');
%  end
%  
%end
%fclose(fid);
%sectionAbbreviation

sectionAbbreviation = [...
    'Drake Pa';...
    'Tasm-Ant';...
    'Afri-Ant';...
    'Antilles';...
    'Mona Pas';...
    'Wind Pas';...
    'FL-Cuba ';...
    'FL-Baham';...
    'Ind Thru';...
    'Agulhas ';...
    'Mozam Ch';...
    'Bering  ';...
];

subplot(2,1,1)
plot(uTransportSection)
axis tight
grid on
xlabel('time counter')
ylabel('transport, Sv')
title('total flow')
for j=1:nSections
  legendText(j) = {[ sectionAbbreviation(j,:) ':' ...
		    num2str(mean(uTransportSection(:,j)),'%5.2f') ]};
end
legend(legendText,'location','NorthEast');

uTransportMinusMean = zeros(size(uTransportSection));
for j=1:nSections
  uTransportMinusMean(:,j) = uTransportSection(:,j) - mean(uTransportSection(:,j));
end
subplot(2,1,2)
plot(uTransportMinusMean)
axis tight
grid on
xlabel('time counter')
ylabel('transport, Sv')
title('flow variability (total-mean)')

   set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	   'PaperPosition',[0.25 0.25 8 8])

   subplot('position',[0 .95 1 .05]); axis off
   text(.005,.7,[ date ]);
