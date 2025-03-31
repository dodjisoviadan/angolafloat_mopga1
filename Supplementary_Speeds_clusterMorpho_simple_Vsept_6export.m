
cd('C:\MOPGA 2022-2023 Dodji\ARGO float Angola')
% Make a choice
msg = "Choose your favorite plot";
opts = ["speed" "attenuation" "AvarS"];
option= menu(msg,opts);
choice=opts(option);
% literature
 xr=[-2.4495 1.4322];          yr=[0.5606 3.0227];
 xn=[-2.3458 1.3198];          yn=[1.0454 2.7045];
xpn=[-1.7752 1.0518];          ypn=[0.0454 3.3636];

XR=[-2.4495 1.4322];          
YR=0.62*XR+log10(132)-1; YR=[-0.3981 2.0085];%[1.6019    4.0085];

MDS_Environment=readtable('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/MDS_Environment.csv');
P1=131:180; P2=181:235; P3=251:310; P4=321:380; P5=416:480; P6=481:545; 
tp_list ={'2021-08-19 to 2021-09-25','2021-10-01 to 2021-11-13','2021-11-13 to 2021-12-15', '2022-01-20 to 2022-03-06', '2022-03-10 to 2022-04-16'};
tp_list ={'2021-07-20 to 2021-08-15','2021-08-19 to 2021-09-20','2021-10-15 to 2021-11-10', '2021-11-20 to 2021-12-15', '2022-01-30 to 2022-02-25', '2022-03-15 to 2022-04-10'};
tp_list ={'2021-07-17 to 2021-08-19','2021-08-22 to 2021-09-24','2021-10-08 to 2021-11-13', '2021-11-23 to 2021-12-31', '2022-01-28 to 2022-03-09', '2022-03-13 to 2022-04-12'};


cd('C:\MOPGA 2022-2023 Dodji\ARGO float Angola\EcotaxaPart_BCG_UVP6_Plot\ClusterCenterQselection\SPEED\AB_or_BV')
A1=dir('*.csv');filenames1={A1.name}; folder1={A1.folder};
filenames1=string([char(folder1'), repmat('\',length(folder1),1), char(filenames1')]);
cd('C:\MOPGA 2022-2023 Dodji\ARGO float Angola\EcotaxaPart_BCG_UVP6_Plot\ClusterCenterQselection\SPEED')
A2=dir('*.csv');  filenames2={A2.name}; folder2={A2.folder};
filenames2=string([char(folder2'), repmat('\',length(folder2),1), char(filenames2')]);
filenames2=filenames2(contains(filenames2,'Q100'));

filenames=[filenames1;filenames2];
filenames=filenames(contains(filenames,'speed')); % consider only speed and attenuation files and remove biovolume and detritus files
filenames=filenames(~contains(filenames,'ParticleBV'));
filenames=filenames(~contains(filenames,'detritus'));
ill=contains(filenames,tp_list);
filenames=filenames(ill);
[n,m]=size(filenames);
[m,n]=size(filenames);

filenames=[filenames(~contains(filenames,'ParticleAB'));filenames(contains(filenames,'ParticleAB'))]; % ramener les particules sans en bas



N=numel(tp_list);
pas=1:N:m;
%fl=filenames(pas)'
period={'d','s','o','^','p','x'};%,'+'
%clus={'r','g','k','[0.6 0 0.2]','m'};
clus5={'b','r','k','g','m'};
clus5={'#DE64B1','#57A7B3',"#F5C710",'k',"#7E2F8E"};
period=repmat(period,1,7)';
clus=repmat(clus5,N,1);
clus=clus(:);
%[rgb]{1 .5 .2}
str = {
    '\color{red}CLUSTER 1', ...
    '\color{green}CLUSTER 2', ...
    '\color{black}CLUSTER 3', ...
    '\color[rgb]{0.6 0 0.2}CLUSTER 4', ...
    '\color{magenta}PARTICLE AB'};
str = {
    '\color{blue}CLUSTER 1', ...
    '\color{red}CLUSTER 2', ...
    '\color{black}CLUSTER 3', ...
    '\color{green}CLUSTER 4', ...
    '\color{magenta}PARTICLE AB'};

str = {
    '\color[rgb]{0.8706,0.3922,0.6941}CLUSTER 1', ...
    '\color[rgb]{0.3412,0.6549,0.7020}CLUSTER 2', ...
    '\color[rgb]{0.9608,0.7804,0.0627}CLUSTER 3', ...
    '\color{black}CLUSTER 4', ...
    '\color[rgb]{0.6350 0.0780 0.1840}PARTICLE AB'};

figure

for sp=1:(N+1) % manage the loop
if sp==(N+1)
    iT=1:m;
    il=contains(filenames,tp_list);
    sp=[sp sp+1];
else
    iT=sp:N:m;
    il=contains(filenames,tp_list(iT(1)));
end

subplot(2,round(N/2+1/2),sp)

SLOPE_P=[];
BINS_P=[];
for it=iT 
if il(it)==1
S=readtable(char(filenames(it)),'Filetype','text','ReadVariableNames',1); % read file and extract informations
slopes(it).S=-S.slope;
bins(it).b=(S.bin_min+S.bin_max)/2;
slope=-S.slope;
attenuation=-S.attenuation;
bin=(S.bin_min+S.bin_max)/2;
sa=slope./attenuation;

bin(isempty(bin))=NaN;
slope(isempty(slope))=NaN;
attenuation(isempty(attenuation))=NaN;
S.bin_min(isempty(S.bin_min))=NaN;
S.bin_max(isempty(S.bin_max))=NaN;

slope(slope<=0)=NaN;
attenuation(attenuation<0)=NaN;
attenuation(isnan(slope))=NaN;
binsa(it).sa=[bin slope attenuation S.bin_min S.bin_max];

SLOPE_P=[SLOPE_P; slope];
BINS_P=[BINS_P; bin];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(choice,'speed')==1  % create plot
  if strcmp(char(clus{it}),clus5{5})==1
    plot(log10(bin),log10(slope),char(period(it)),'Color', char(clus{it}),'MarkerSize',12,'MarkerFaceColor',char(clus{it}), 'LineWidth',1)
  else
    plot(log10(bin),log10(slope),char(period(it)),'Color', char(clus{it}),'MarkerSize',12,'LineWidth',1 )
  end
hold on
elseif strcmp(choice,'attenuation')==1
  if strcmp(char(clus{it}),clus5{5})==1
    plot(log10(bin),(attenuation),char(period(it)),'Color', char(clus{it}),'MarkerSize',12,'MarkerFaceColor',char(clus{it}), 'LineWidth',1)
  else
    plot(log10(bin),(attenuation),char(period(it)),'Color', char(clus{it}),'MarkerSize',12, 'LineWidth',1 )
  end
  hold on
else
  if strcmp(char(clus{it}),clus5{5})==1
    plot(log10(slope),(attenuation),char(period(it)),'Color', char(clus{it}),'MarkerSize',12,'MarkerFaceColor',char(clus{it}), 'LineWidth',1)
  else
    plot(log10(slope),(attenuation),char(period(it)),'Color', char(clus{it}),'MarkerSize',12, 'LineWidth',1 )
  end  
hold on
end

end
end

if strcmp(choice,'speed')==1     % create plot and paramater
plot(xr,yr,'-','color',[0.9290, 0.6940, 0.1250])
plot(xn,yn,'-k')
plot(xpn,ypn,'--k')
plot(XR,YR,'--','color','g','LineWidth',1.5)

xlabel 'Log10(Size Classes,mm)'
ylabel 'Log10(Speed m/d)'
xlim([0 3.5])
ylim([0 90])
xlim([-2.5 1.5])
xlim([-1.1 0.55])
ylim([0 3.5])
title 'Speeds variation with size classes and Detritus or Clusters'
title 'Speeds variation with size'
title ''
SLOPES_positive(sp(1)).S=SLOPE_P;
elseif strcmp(choice,'attenuation')== 1                  % create plot paramete
xlabel 'Log10(Size Classes,mm)'
ylabel 'attenuation abundance'
xlim([0 3.5])
xlim([-2.5 1.5])
xlim([-1.1 0.55])
ylim([0 2.5])
title 'Attenation variation with size classes and Detritus or Clusters'
title 'Attenation variation with size'
title ''
else
   xlabel 'Log10(Speed m/d)'
ylabel 'attenuation abundance'
xlim([0 3.5])
ylim([0 2.5])
title 'Attenuation variation with Speeds'
title ''

end

if sp(1)==(N+1)
    lg=legend(tp_list,'Location','northwest');
else
    lg=legend(tp_list(sp),'Location','northeast');
    if sp==1 & strcmp(choice,'speed')~=1
        % add text
         t = text(log10(0.5),1,str,'Position',[-1.94056783940832 1.83474727434081 0], 'fontweight', 'bold','Position',[7.28392880488698 -1.55234949985274 0] ); 
    elseif sp==1 & strcmp(choice,'speed')==1
        % add text
t = text(log10(0.5),log10(80),str,'Position',[-1.94056783940832 1.83474727434081 0], 'fontweight', 'bold','Position',[7.28392880488698 -1.55234949985274 0] );
text(log10(0.5),log10(80),'OLS','Position',[-1.65885348540287 -0.572365825726559 0]);
text(log10(0.5),log10(80),'Hierarchical','Position',[-1.64879792449464 -0.85100900844993 0]);
text(log10(0.5),log10(80),'BETA=1.17','Position',[-1.65502845408966 -0.358736804175026 0]);
text(log10(0.5),log10(80),' From Cael et al.2020','Position',[-2.00171022554639 -0.0937237493447391 0]);
text(log10(0.5),log10(80),' From Kriest et al.2022','Position',[-2.00171022554639 -0.0937237493447391 0]);

% Create line
annotation('line',[0.0188216039279869 0.0510416666666667],...
    [0.520889043053282 0.520926090828138]);

% Create line
annotation('line',[0.0189707787779596 0.0515625],...
    [0.542908550601585 0.542297417631345],'LineStyle','--');

% Create line
annotation('line',[0.0177083333333333 0.0505131614839062],...
    [0.496883348174533 0.496726997867009],...
    'Color',[1 0.411764705882353 0.16078431372549]);
% Create line
annotation('line',[0.0177083333333333 0.0505131614839062],...
    [0.496883348174533 0.496726997867009],'LineStyle','--',...
    'Color',[0 1 0]);


    end
end
%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontSize',12)
end
%0/aaxa
cd('C:\MOPGA 2022-2023 Dodji\ARGO float Angola\EcotaxaPart_BCG_UVP6_Plot\ClusterCenterQselection\SPEED\AOUT2024')

f=gcf;
if strcmp(choice,'speed')==1  
    exportgraphics(f,'Fev2024Speed.png','Resolution',1200)
    exportgraphics(f,'Fev2024Speed.pdf','Resolution',1200)
elseif strcmp(choice,'attenuation')==1  
    exportgraphics(f,'Fev2024Attenuation.png','Resolution',1200)
    exportgraphics(f,'Fev2024Attenuation.pdf','Resolution',1200)
else
    exportgraphics(f,'Fev2024Attenuation_vs_Speed.png','Resolution',1200)
    exportgraphics(f,'Fev2024Attenuation_vs_Speed.pdf','Resolution',1200)
end

%% SPEED AND ATTENUATION ON ALL COLUMNS OF ONE CLUSTER
dk=1:6:24;
binsak_cluters=[];

for k=1:4
    binsak=[];
    for ik=dk(k):dk(k)+5
        binsak=[binsak;binsa(ik).sa];
    end
    binsak=binsak(0.7290<=binsak(:,1)<=1.46,:);
    binsak=binsak(binsak(:,2)>0,:);

    cz=0;
    for sz=[0.7290 0.9165 1.1550 1.46]
        cz=cz+1;
        binsaksize(cz).size=binsak(binsak(:,1)==sz,:);
        bslp=binsak(binsak(:,1)==sz,2);
        bslp=bslp(~isoutlier(bslp,'median'));
        batt=binsak(binsak(:,1)==sz,3);
        batt=batt(~isoutlier(batt,'median'));
        binsaksize(cz).sizemean=[mean(bslp) mean(batt)];
        binsaksize(cz).sizeSTD=[std(bslp) std(batt)];
        binsaksize(cz).sizeperiod=[(binsak(binsak(:,1)==sz,2)) (binsak(binsak(:,1)==sz,3))];
    end

binsak_cluters(k).cluster=binsaksize;
end
0/0
