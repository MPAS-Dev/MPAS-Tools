% Eulerian velocity from prognostic momentum equation
%var_name = {'avgNormalVelocity'};
% total transport velocity
var_name = {'avgNormalTransportVelocity'}

% plot transport results for gm

sectionText = {
'Drake Passage, S Ocean -56  to -63 lat,  68W lon, section A21,  140+/- 6 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Tasmania-Ant, S Ocean  -44  to -66 lat, 140E lon, section P12,  157+/-10 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Africa-Ant, S Ocean    -31.3to -70 lat,  30E lon, section  I6,           Sv in Ganachaud99thesis                        ',...
'Antilles Inflow, Carib.                                       -18.4+/-4.7Sv in Johns_ea02dsr                            '...
'Mona Passage, Caribbian                                        -2.6+/-1.2Sv in Johns_ea02dsr                            '...
'Windward Passage, Carib                                        -7.0      Sv in Nelepo_ea76sr, Roemmich81jgr             '...
'Florida-Cuba, Caribbian                                        31.5+/-1.5Sv in Johns_ea02dsr, 32.3+/-3.2Sv Larsen92rslpt'...
'Florida-Bahamas, Carib. 27 lat, -80  to -78.8lon,              31.5+/-1.5Sv in Johns_ea02dsr, 32.3+/-3.2Sv Larsen92rslpt'...
'Indonesian Throughflow, -9  to -18 lat, 116E lon, section J89,  -16+/- 5 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Agulhas                                                         -70+/-20 Sv in Bryden_Beal01dsr                         ',...
'Mozambique Channel,    -25 lat,  35  to  44E lon, section I4 ,  -14+/- 6 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Bering Strait, Arctic                                          0.83+/-0.66Sv in Roach_ea95jgr                           '...
    };

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

nSims = 6;
ntrans = 12;

dir = 'm91';
abc = 'klmnop';
kappa = [80 100 200 500 1000 2000];

data_mean = zeros(nSims,ntrans);
data_std = zeros(size(data_mean));
data_max = zeros(size(data_mean));
data_min = zeros(size(data_mean));


for iSim=1:nSims
  filename = ['data/' dir abc(iSim) '_' char(var_name) '_small_data_file.mat'];
  load(filename,'mean_transport','std_transport','min_transport','max_transport');

  data_mean(iSim,:) =  mean_transport;
  data_std(iSim,:) =  std_transport;
  data_max(iSim,:) =  max_transport;
  data_min(iSim,:) =  min_transport;
end


for iTrans = 1:3
  figure(iTrans); clf
  semilogx(kappa(2:nSims),data_mean(2:nSims,iTrans),'*-k')
  hold on
  semilogx(kappa(2:nSims),data_mean(2:nSims,iTrans)+data_std(2:nSims,iTrans),'*-r')
  semilogx(kappa(2:nSims),data_min(2:nSims,iTrans),'*-g')
  semilogx(kappa(2:nSims),data_mean(2:nSims,iTrans)-data_std(2:nSims,iTrans),'*-r')
  semilogx(kappa(2:nSims),data_max(2:nSims,iTrans),'*-g')

  h=semilogx(kappa(1),data_mean(1,iTrans),'s-k');
  set(h,'MarkerEdgeColor','none','MarkerFaceColor','k')
  h=semilogx(kappa(1),data_mean(1,iTrans)+data_std(1,iTrans),'s-r');
  set(h,'MarkerEdgeColor','none','MarkerFaceColor','r')
  h=semilogx(kappa(1),data_min(1,iTrans),'s-g');
  set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
  h=semilogx(kappa(1),data_mean(1,iTrans)-data_std(1,iTrans),'s-r');
  set(h,'MarkerEdgeColor','none','MarkerFaceColor','r')
  h=semilogx(kappa(1),data_max(1,iTrans),'s-g');
  set(h,'MarkerEdgeColor','none','MarkerFaceColor','g')
  
  grid on
  hold on
  ttxt = char(sectionText(iTrans))
  title(ttxt(1:23))
  xlabel('GM kappa value, m^2/s.  Far left is GM off')
  ylabel('transport, Sv, 3-year mean')
  legend('3 year mean','mean+/-monthly std dev','min/max of monthly means','Location','SouthWest' )
  xlim([70 2500])

  % for log-log plot comparison to Dana95
  %axis([1e2 1e4 50 200])
  %set(gca,'YTick',[50 100 200])
  
  set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4 3])
  fig=[char(dir) abc(1:nSims) '_'  char(var_name) '_' num2str(iTrans)];
  print('-depsc2',['f/' fig '.eps']);
  print('-djpeg',['f/' fig '.jpg']);
  unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

end
