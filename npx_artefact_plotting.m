
%this is the code for producing the analysis and figures for the manuscript
%"Overcoming off-target optical stimulation-evoked cortical activity in the
%mouse brain in vivo"; Simon Weiler, Mateo Velez-Fort  and Troy W. Margrie

%this uses only data obtained in the year 2024 
%% DEPENDENCIES, PLEASE READ

%Please get the follwing functions: 
%1) uipickfiles: https://uk.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids
%2) shadedErrorBar: https://uk.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
%3) violinplot: https://github.com/bastibe/Violinplot-Matlab

%4) npx_example_plotter.m: same github path

%needs Data Structure : final_structure.mat which is on the github path

%ID for cortical layers in VISp (primary visual cortex)
%593= VISp1;821= VISp2/3; 721= VISp4; 778= VISp5; 33= VISp6a; 305= VISp6b;

%% Sessions
%1) blue darkness (1-5 laser intensities)
%2) orange darkness (1-5 laser intensities)
%3) red darkness (1-5 laser intensities)

%4) red ambient light 20 lux (1-5 laser intensities)
%5) red ambient light 40 lux (1-5 laser intensities)
%6) red ambient light 80 lux (1-5 laser intensities)

%7) blue ambient light 20 lux (1-5 laser intensities)
%8) blue ambient light 40 lux (1-5 laser intensities)
%9) blue ambient light 80 lux (1-5 laser intensities)

%10) orange ambient light 20 lux (1-5 laser intensities)
%11) orange ambient light 40 lux (1-5 laser intensities)
%12) orange ambient light 80 lux (1-5 laser intensities)
%% Load data structure  
clear all;clc;
%using uipickfiles
str_neuropixels    = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\npx_artefact\output';%change accordingly
folder_list = uipickfiles('FilterSpec',str_neuropixels);load(char(folder_list));
%folder location where to save figure panels, change accordingly
save_folder='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\npx_artefact\figure_panels';
%% Set parameters
%arial as default font
set(0, 'DefaultAxesFontName', 'Arial'); 
%colors for laser lights 
blue_c=([0 191 255]/256);orange_c=([255 165 0]/256);red_c=[1 0 0];
%time before and after laser stimulus in seconds
t_pre=0.2;t_post=0.8;
%laser powers 
l_power = [1 2.5 5 10 15];
%animal lenghts 
for i=1:size(neuropixel_data,2)
    units_pa(i)=length(neuropixel_data(i).spike_type);
end
animal_s=cumsum(units_pa);a_units={1:animal_s(1); animal_s(1)+1:animal_s(2); animal_s(2)+1:animal_s(3); animal_s(3)+1:animal_s(4)};
animal_idx=[ones(units_pa(1),1) ;2*ones(units_pa(2),1)...
    ;3*ones(units_pa(3),1) ;4*ones(units_pa(4),1)];
l23id=821;
l4id=721;
l5id=778;
l6aid=33;
l6bid=305;
%% read out zeta values and delta HZ values and SDF traces for different laser intensities and color for a given session across all animals (concatenate)
%BLUE (see above for which sessions are blue) 
sess_i=[1 7 8 9];
sdfs_blue=[];sdfm_blue=[];MIs_blue=[];MIm_blue=[];zetas_blue=[];zetam_blue=[];zetasn_blue=[];
for k=1:length(sess_i)
    mods_idx=[];modm_idx=[];zetas=[];zetam=[];zeta_noise=[];
    for l_inten=1:5;
    [temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8 temp9 temp10 temp11] = across_animals(neuropixel_data,sess_i(k),l_inten);
    sdfs_blue(:,:,l_inten,k)=[temp1];sdfm_blue(:,:,l_inten,k)=[temp2];
    mods_idx=[mods_idx temp3(:,2)];modm_idx=[modm_idx temp4(:,2)];
    zetas=[zetas temp5];zetam=[zetam temp6];zeta_noise=[zeta_noise temp7];
    end
    MIs_blue(:,:,k)=mods_idx;MIm_blue(:,:,k)=modm_idx;
    zetas_blue(:,:,k)=zetas;zetam_blue(:,:,k)=zetam;zetasn_blue(:,:,k)=zeta_noise;
end


%ORANGE (see above for which sessions are orange) 
sess_i=[2 10 11 12];
sdfs_orange=[];sdfm_orange=[];MIs_orange=[];MIm_orange=[];zetas_orange=[];zetam_orange=[];zetasn_orange=[];
for k=1:length(sess_i)
    mods_idx=[];modm_idx=[];zetas=[];zetam=[];zeta_noise=[];
    for l_inten=1:5;
    [temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8 temp9 temp10 temp11] = across_animals(neuropixel_data,sess_i(k),l_inten);
    sdfs_orange(:,:,l_inten,k)=[temp1];sdfm_orange(:,:,l_inten,k)=[temp2];
    mods_idx=[mods_idx temp3(:,2)];modm_idx=[modm_idx temp4(:,2)];
    zetas=[zetas temp5];zetam=[zetam temp6];zeta_noise=[zeta_noise temp7];
    end
    MIs_orange(:,:,k)=mods_idx;MIm_orange(:,:,k)=modm_idx;
    zetas_orange(:,:,k)=zetas;zetam_orange(:,:,k)=zetam;zetasn_orange(:,:,k)=zeta_noise;
end

%RED (see above for which sessions are red) 
sess_i=[3 4 5 6];
sdfs_red=[];sdfm_red=[];MIs_red=[];MIm_red=[];zetas_red=[];zetam_red=[];zetasn_red=[];
for k=1:length(sess_i)
    mods_idx=[];modm_idx=[];zetas=[];zetam=[];zeta_noise=[];
    for l_inten=1:5;
    [temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8 temp9 temp10 temp11] = across_animals(neuropixel_data,sess_i(k),l_inten);
    sdfs_red(:,:,l_inten,k)=[temp1];sdfm_red(:,:,l_inten,k)=[temp2];
    mods_idx=[mods_idx temp3(:,2)];modm_idx=[modm_idx temp4(:,2)];
    zetas=[zetas temp5];zetam=[zetam temp6];zeta_noise=[zeta_noise temp7];
    end
    MIs_red(:,:,k)=mods_idx;MIm_red(:,:,k)=modm_idx;
    zetas_red(:,:,k)=zetas;zetam_red(:,:,k)=zetam;zetasn_red(:,:,k)=zeta_noise;
end
%depth ID and values 
depth_id=temp8;
depth_nr=temp9;
depth_idm=temp10;
depth_nrm=temp11;
%% Following plots the example traces in darkness and ambient light with different laser intensities for
%FIG 1B and FIG 2B

%% Blue lowest power 
unit_id=104;sess_i=[];inten_i=[];sess_i=1;sess_r=1;inten_i=1;raster_ind=1:10;color_s=blue_c;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];

plot_pre=0.2;plot_post=0.8;
color_s=blue_c;
average_sdf=sdfs_blue(:,:,inten_i,sess_i);
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'1 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_blue_lowest_power.pdf');
%% Blue highest power 
unit_id=104;sess_i=[];inten_i=[];sess_i=1;sess_r=1;inten_i=5;raster_ind=41:50;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_blue(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=blue_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'15 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_blue_highest_power.pdf');
%% Orange lowest power 
unit_id=705;sess_i=[];inten_i=[];sess_i=1;sess_r=2;inten_i=1;raster_ind=1:10;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_orange(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=orange_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'1 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_orange_lowest_power.pdf');
%% Orange highest power 
unit_id=705;sess_i=[];inten_i=[];sess_i=1;sess_r=2;inten_i=5;raster_ind=41:50;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_orange(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=orange_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'15 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_orange_highest_power.pdf');
%% Red lowest power 
unit_id=651;sess_i=[];inten_i=[];sess_i=1;sess_r=3;inten_i=1;raster_ind=1:10;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_red(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=red_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'1 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_red_lowest_power.pdf');
%% Red highest power 
unit_id=651;sess_i=[];inten_i=[];sess_i=1;sess_r=3;inten_i=5;raster_ind=41:50;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_red(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=red_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'1 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_red_highest_power.pdf');
%% Same units AMBIENT lights on
%% Blue 20 lux light on max power 
unit_id=104;sess_i=[];inten_i=[];sess_i=2;sess_r=7;inten_i=5;raster_ind=41:50;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_blue(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=blue_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'15 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_blue_highest_power_AMBIENT.pdf');
%% Orange 20 lux light on highest power 
unit_id=705;sess_i=[];inten_i=[];sess_i=2;sess_r=10;inten_i=5;raster_ind=41:50;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_orange(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=orange_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'15 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_orange_highest_power_AMBIENT.pdf');
%% Red 20 lux light on highest power 
unit_id=651;sess_i=[];inten_i=[];sess_i=2;sess_r=4;inten_i=5;raster_ind=41:50;
raster=[];average_sdf=[];
raster=[neuropixel_data(1).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(2).session_aligned(sess_r).raster(raster_ind,:)...
    neuropixel_data(3).session_aligned(sess_r).raster(raster_ind,:) neuropixel_data(4).session_aligned(sess_r).raster(raster_ind,:)];
average_sdf=sdfs_red(:,:,inten_i,sess_i);

plot_pre=0.2;plot_post=0.8;
color_s=red_c;
npx_example_plotter(unit_id,raster,average_sdf,plot_pre,plot_post,color_s);
text(-0.05,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,'15 mW','FontWeight','normal','Color',color_s);
text(-0.35,max(get(gca, 'ylim'))+max(get(gca, 'ylim'))/2,['unit ' num2str(unit_id)],'FontWeight','normal');
% save
cd(save_folder);saveas(gcf, 'example_unit_red_highest_power_AMBIENT.pdf');
%% Plot on top of each other all units , different intesities for red, orange and blue 
%% Darkness 15 mW Figure 1C top
% Plot SDF traces over time 
lux_m='0';sess_i=1;inten_i=5;scalef=1;
temp1=[];temp2=[];temp3=[];
%only single units 
temp1=sdfs_blue(:,:,inten_i,sess_i);
temp2=sdfs_orange(:,:,inten_i,sess_i);
temp3=sdfs_red(:,:,inten_i,sess_i);
stim_on=20;

average_errorbar_sdf(temp1,temp2,temp3,sess_i,inten_i,scalef,l_power,stim_on,blue_c,orange_c,red_c,lux_m)

% save
cd(save_folder);saveas(gcf, 'average_dark_15mW.pdf');
%% Darkness 1 mW Figure 1C bottom
% Plot SDF traces over time 
lux_m='0';sess_i=1;inten_i=1;scalef=1;
temp1=[];temp2=[];temp3=[];
%only single units 
temp1=sdfs_blue(:,:,inten_i,sess_i);
temp2=sdfs_orange(:,:,inten_i,sess_i);
temp3=sdfs_red(:,:,inten_i,sess_i);
stim_on=20;

average_errorbar_sdf(temp1,temp2,temp3,sess_i,inten_i,scalef,l_power,stim_on,blue_c,orange_c,red_c,lux_m)

% save
cd(save_folder);saveas(gcf, 'average_dark_1mW.pdf');
%% Ambient light 20 lux,  15 mW Figure 1E
% Plot SDF traces over time 
lux_m='20';sess_i=2;inten_i=5;scalef=1;
temp1=[];temp2=[];temp3=[];
%only single units 
temp1=sdfs_blue(:,:,inten_i,sess_i);
temp2=sdfs_orange(:,:,inten_i,sess_i);
temp3=sdfs_red(:,:,inten_i,sess_i);
stim_on=20;

average_errorbar_sdf(temp1,temp2,temp3,sess_i,inten_i,scalef,l_power,stim_on,blue_c,orange_c,red_c,lux_m)

% save
cd(save_folder);saveas(gcf, 'average_20lux_15mW.pdf');
%%  Ambient light 20 lux,  2.5 mW Figure 1E
% Plot SDF traces over time 
lux_m='20';sess_i=2;inten_i=2;scalef=1;
temp1=[];temp2=[];temp3=[];
%only single units 
temp1=sdfs_blue(:,:,inten_i,sess_i);
temp2=sdfs_orange(:,:,inten_i,sess_i);
temp3=sdfs_red(:,:,inten_i,sess_i);
stim_on=20;

average_errorbar_sdf(temp1,temp2,temp3,sess_i,inten_i,scalef,l_power,stim_on,blue_c,orange_c,red_c,lux_m)

% save
cd(save_folder);saveas(gcf, 'average_20lux_2.5mW.pdf');
%% Baseline vs light response PERMUTATION test
%change session and intesity accordingly 
sess_i=1;intens_i=1;
baser=1:200;
respr=201:700;
temp_sdf={sdfs_blue,sdfs_orange,sdfs_red};
tmp_zeta={zetas_blue,zetas_orange,zetas_red};
p_all=[];p=[];
for k=1:10
for i=1:3
    sdf=[];sdf=temp_sdf{i};
    zta=[];zta=tmp_zeta{i};trz=[];tr1=[];tr2=[];
    tr1=nanmean(sdf(:,baser,intens_i,sess_i),2);
    tr2=nanmean(sdf(:,respr,intens_i,sess_i),2);
    trz=zta(:,intens_i,sess_i);
 %[r(i) p(i)]=signrank(tr1(~isnan(trz)), tr2(~isnan(trz)));
 %[p(i), observeddifference, effectsize] = permutationTest(tr1(~isnan(trz)), tr2(~isnan(trz)), 1000);
 [p(i), observeddifference, effectsize] = permutationTest(tr1, tr2, 1000);

mb(i)=nanmean(tr1(~isnan(trz)));
mr(i)=nanmean(tr2(~isnan(trz)));
mb_sem(i)=nanstd(tr1(~isnan(trz)))/sqrt(length(tr1(~isnan(trz))));
mr_sem(i)=nanstd(tr2(~isnan(trz)))/sqrt(length(tr2(~isnan(trz))));

    %[h(i) p] = adtest(nanmean(sdf(:,respr,intens_i,sess_i),2));
    %[r(i) p(i)]=signrank(nanmean(sdf(:,baser,intens_i,sess_i),2), nanmean(sdf(:,respr,intens_i,sess_i),2));
  % [r(i) p(i)]=ranksum(allbase(:), nanmean(sdf(:,respr,intens_i,sess_i),2));
end
p_all(k,:)=p;
end
mb
mb_sem
mr
mr_sem

%% Average percentage across laser intensities darkness only, Bar plots in Fig 1B
color_all=[];
color_all=[blue_c; orange_c; red_c]
sess_i=1;
%single or multi
zeta_temp1=[];zeta_temp2=[];zeta_temp3=[];temp1a=[];temp2a=[];temp3a=[];
zeta_temp1=zetas_blue;zeta_temp2=zetas_orange;zeta_temp3=zetas_red;
data_m=[];temp1=[];temp2=[];temp3=[];
temp1=(sum(zeta_temp1(:,:,sess_i)<0.05)/length(zeta_temp1))*100
temp2=(sum(zeta_temp2(:,:,sess_i)<0.05)/length(zeta_temp2))*100
temp3=(sum(zeta_temp3(:,:,sess_i)<0.05)/length(zeta_temp3))*100
temp1a=(sum(zeta_temp1(:,:,sess_i)<0.05))
temp2a=(sum(zeta_temp2(:,:,sess_i)<0.05))
temp3a=(sum(zeta_temp3(:,:,sess_i)<0.05))

fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [900, 400 ,900, 300]);t=tiledlayout("horizontal",'TileSpacing','Compact');
temp_all=[temp1 ;temp2 ;temp3];
for i=1:3
nexttile
b1=bar(1,nanmean(temp_all(i,:)),0.2);hold on;ylim([0 50]);
box off;set(gca,'TickDir','out');ylabel('Modulated units (%)');xlabel('Fiber laser power (mW)');
set(b1,'ShowBaseLine','off');b1(1).FaceColor=color_all(i,:);
h = gca;h.XAxis.Visible = 'off';set(gca,'FontSize',12);
sc1=scatter(ones(1,length(temp_all(i,:))),temp_all(i,:),45);hold on;sc1.MarkerEdgeColor=[0.5 0.5 0.5];sc1.MarkerEdgeAlpha=0.8;%sc1.MarkerFaceColor=[0.5 0.5 0.5];
%sc1.MarkerFaceAlpha=0.8;
errorbar(1,nanmean(temp_all(i,:)),nanstd(temp_all(i,:))/sqrt(length(temp_all(i,:))),'Color','k','LineWidth',1.2);

end
%set(b1,'ShowBaseLine','off');b1(1).FaceColor=blue_c;b1(2).FaceColor=orange_c;b1(3).FaceColor=red_c;
cd(save_folder);saveas(gcf, 'percentage_dark_average_across_intensities.pdf');
%% Average percentage across laser intensities darkness + LIGHT, bar plots in Fig 2B
color_all=[];
color_all=[blue_c; orange_c; red_c]

%single or multi
zeta_temp1=[];zeta_temp2=[];zeta_temp3=[];temp1a=[];temp2a=[];temp3a=[];
zeta_temp1=zetas_blue;zeta_temp2=zetas_orange;zeta_temp3=zetas_red;
data_m=[];temp1=[];temp2=[];temp3=[];temp4=[];temp5=[];temp6=[];

temp1=(sum(zeta_temp1(:,:,1)<0.05)/length(zeta_temp1))*100
temp2=(sum(zeta_temp2(:,:,1)<0.05)/length(zeta_temp2))*100
temp3=(sum(zeta_temp3(:,:,1)<0.05)/length(zeta_temp3))*100

temp4=(sum(zeta_temp1(:,:,2)<0.05)/length(zeta_temp1))*100
temp5=(sum(zeta_temp2(:,:,2)<0.05)/length(zeta_temp2))*100
temp6=(sum(zeta_temp3(:,:,2)<0.05)/length(zeta_temp3))*100

fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [900, 400 ,900, 300]);t=tiledlayout("horizontal",'TileSpacing','Compact');
temp_all=[temp1 ;temp2 ;temp3];
temp_all2=[temp4 ;temp5 ;temp6];
for i=1:3
nexttile
b1=bar(1,nanmean(temp_all(i,:)),0.2);hold on;ylim([0 50]);
b2=bar(1.3,nanmean(temp_all2(i,:)),0.2);hold on;ylim([0 50]);
box off;set(gca,'TickDir','out');ylabel('Modulated units (%)');xlabel('Fiber laser power (mW)');
set(b1,'ShowBaseLine','off');b1(1).FaceColor=color_all(i,:);
set(b2,'ShowBaseLine','off');b2(1).FaceColor=color_all(i,:);
h = gca;h.XAxis.Visible = 'off';set(gca,'FontSize',12);
sc1=scatter(ones(1,length(temp_all(i,:))),temp_all(i,:),45);hold on;sc1.MarkerEdgeColor=[0.5 0.5 0.5];sc1.MarkerEdgeAlpha=0.8;%sc1.MarkerFaceColor=[0.5 0.5 0.5];
%sc1.MarkerFaceAlpha=0.8;
sc2=scatter(ones(1,length(temp_all2(i,:)))*1.3,temp_all2(i,:),45);hold on;sc2.MarkerEdgeColor=[0.5 0.5 0.5];sc2.MarkerEdgeAlpha=0.8;%sc2.MarkerFaceColor=[0.5 0.5 0.5];
%sc2.MarkerFaceAlpha=0.8;
errorbar(1,nanmean(temp_all(i,:)),nanstd(temp_all(i,:))/sqrt(length(temp_all(i,:))),'Color','k','LineWidth',1.2);
errorbar(1.3,nanmean(temp_all2(i,:)),nanstd(temp_all2(i,:))/sqrt(length(temp_all2(i,:))),'Color','k','LineWidth',1.2);
end
%set(b1,'ShowBaseLine','off');b1(1).FaceColor=blue_c;b1(2).FaceColor=orange_c;b1(3).FaceColor=red_c;
cd(save_folder);saveas(gcf, 'percentage_dark_averageVSLIGHT_across_intensities.pdf');
%% compare darkness and ambient light fraction 
zeta_temp1=[];
zeta_temp1=zetas_red;
temp1=[];temp2=[];
sess_ia=1;sess_ib=2;
laser_i=5;
temp1=(sum(zeta_temp1(:,laser_i,sess_ia)<0.05))
temp2=(sum(zeta_temp1(:,laser_i,sess_ib)<0.05))
 x1=[];x1=[temp1,length(MIs_blue)-temp1;temp2,length(MIs_blue)-temp2];
 [h,p,stats] = fishertest(x1);
%% fisher test all
for i=1:5
    x1=[];x1=[temp1a(i),length(MIs_blue)-temp1a(i);temp3a(i),length(MIs_blue)-temp3a(i)];
    x2=[];x2=[temp2a(i),length(MIs_blue)-temp2a(i);temp3a(i),length(MIs_blue)-temp3a(i)];
    [h,p_rb(i),stats] = fishertest(x1);
    [h,p_ro(i),stats] = fishertest(x2);
end
%% average across 5 intensitits
p=[];
[p,tbl,stats] = anova1([temp1' temp2' temp3'],[])
figure;c = multcompare(stats);


%% Delta Hz change (including significant and non signficant ones) -> Everything included 
sess_i=1;
data_m=[];temp1=[];temp2=[];temp3=[];
mi_temp1=[];mi_temp2=[];mi_temp3=[];
mi_temp1=MIs_blue;mi_temp2=MIs_orange;mi_temp3=MIs_red;
% zeta1=[];zeta2=[];zeta3=[];
%  zeta1=zetas_blue(:,intens_i,sess_i);

temp1=nanmean(mi_temp1(:,:,sess_i));
temp1_err=nanstd(mi_temp1(:,:,sess_i))/sqrt(length(mi_temp1));
temp2=nanmean(mi_temp2(:,:,sess_i));
temp2_err=nanstd(mi_temp2(:,:,sess_i))/sqrt(length(mi_temp2));
temp3=nanmean(mi_temp3(:,:,sess_i));
temp3_err=nanstd(mi_temp3(:,:,sess_i))/sqrt(length(mi_temp3));

fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [900, 400 ,350, 300]);
h=errorbar(temp1,temp1_err, '-o','Color', blue_c,'LineWidth', 1,'CapSize',5,"MarkerFaceColor",blue_c);hold on;
h=errorbar(temp2,temp2_err, '-o','Color', orange_c,'LineWidth', 1,'CapSize',5,"MarkerFaceColor",orange_c);hold on;
h=errorbar(temp3,temp3_err, '-o','Color', red_c,'LineWidth', 1,'CapSize',5,"MarkerFaceColor",red_c);
xlim([1 5.5])
ylim([0 3])
box off;set(gca,'TickDir','out');
ylabel('\Delta Firing rate (Hz)');xlabel('Fiber laser power (mW)');
offsetAxes;
xticklabels({'1','2.5','5','10','15'});set(gca,'FontSize',12);
hold on;text(4.25,0.75,'637 nm','Color','r','FontWeight','bold','FontSize',12);
hold on;text(4.25,0.45,'594 nm','Color',([255 165 0]/256),'FontWeight','bold','FontSize',12);
hold on;text(4.25,0.15,'473 nm ','Color',([0 191 255]/256),'FontWeight','bold','FontSize',12);

hold on;text(4.8,3,'*** ','Color','k','FontWeight','bold','FontSize',18);
hold on;text(3.8,3,'*** ','Color','k','FontWeight','bold','FontSize',18);
hold on;text(2.8,3,'*** ','Color','k','FontWeight','bold','FontSize',18);
hold on;text(1.8,3,'*** ','Color','k','FontWeight','bold','FontSize',18);
hold on;text(0.8,3,'*** ','Color','k','FontWeight','bold','FontSize',18);

%hold on;text(4.8,2.75,'* ','Color',[0.5 0.5 0.5],'FontWeight','bold','FontSize',18);
hold on;text(3.8,2.75,'*** ','Color',[0.5 0.5 0.5],'FontWeight','bold','FontSize',18);
%hold on;text(1.8,2.75,'*** ','Color',[0.5 0.5 0.5],'FontWeight','bold','FontSize',18);
hold on;text(0.8,2.75,'*** ','Color',[0.5 0.5 0.5],'FontWeight','bold','FontSize',18);

hold on;line([5.4 5.4],[temp1(5) temp3(5)],'Color','k');
hold on;line([5.35 5.4],[temp1(5) temp1(5)],'Color','k');
hold on;line([5.35 5.4],[temp3(5) temp3(5)],'Color','k');

hold on;line([5.2 5.2],[temp2(5) temp3(5)],'Color',[0.5 0.5 0.5]);
hold on;line([5.15 5.2],[temp2(5) temp2(5)],'Color',[0.5 0.5 0.5]);
hold on;line([5.15 5.2],[temp3(5) temp3(5)],'Color',[0.5 0.5 0.5]);
cd(save_folder);saveas(gcf, 'delta_dark.pdf');
temp1
temp2
temp3
temp1_err
temp2_err
temp3_err
%% Statistcs 
h=[];stats=[];pfried=[];c=[];
i=i;sess_i=1;
tempo=[];tempo=[mi_temp1(:,i,sess_i) mi_temp2(:,i,sess_i) mi_temp3(:,i,sess_i)];
[pfried tbl stats] = friedman(tempo)
figure;c = multcompare(stats);
%[h(i) ] = adtest(tempo(:,1))

% r_b(i)=signrank(tempo(:,1),tempo(:,3));
% r_o(i)=signrank(tempo(:,2),tempo(:,3));
% o_b(i)=signrank(tempo(:,2),tempo(:,1));



%% Time onset BLUE, ORANGE, RED HIGHEST INTENSITY
p1=[];p2=[];p3=[];p1m=[];p2m=[];p3m=[];
temp1=[];temp2=[];temp3=[];temp1m=[];temp2m=[];temp3m=[];
inten_i=5;
sess_i=1;

%select zetas for SUA
p1=zetas_blue(:,inten_i,sess_i);
p2=zetas_orange(:,inten_i,sess_i);
p3=zetas_red(:,inten_i,sess_i);
%select zetas for MUA
p1m=zetam_blue(:,inten_i,sess_i);
p2m=zetam_orange(:,inten_i,sess_i);
p3m=zetam_red(:,inten_i,sess_i);

%only include <0.05 p values 
%SUA
temp1=sdfs_blue(p1<0.05,:,inten_i,sess_i);
temp2=sdfs_orange(p2<0.05,:,inten_i,sess_i);
temp3=sdfs_red(p3<0.05,:,inten_i,sess_i);
%MUA
temp1m=sdfs_blue(p1m<0.05,:,inten_i,sess_i);
temp2m=sdfs_orange(p2m<0.05,:,inten_i,sess_i);
temp3m=sdfs_red(p3m<0.05,:,inten_i,sess_i);
%DEPTH ID BLUE and animal 
blue_depth=depth_id(p1<0.05);
blue_depthm=depth_idm(p1m<0.05);
blue_a=animal_idx(p1<0.05);
blue_am=animal_idx(p1m<0.05);
%DEPTH ID ORANGE animal 
orange_depth=depth_id(p2<0.05);
orange_depthm=depth_idm(p2m<0.05);
orange_a=animal_idx(p2<0.05);
orange_am=animal_idx(p2m<0.05);
%DEPTH ID RED and animal 
red_depth=depth_id(p3<0.05);
red_depthm=depth_idm(p3m<0.05);
red_a=animal_idx(p3<0.05);
red_am=animal_idx(p3m<0.05);


%% ONSET TIMING 
%BLUE
cd('I:\neuropixels_artefact\onsets\blue');
bs=1:200;
act=201:600;
fc=3;
time_onsetb=[];
for i=1:size(temp1,1)
trac=[];
trac=temp1(i,:);
[time_onsetb(i)] = onset_detection(trac,bs,act,fc);
%title(['unit ' num2str(i)]);saveas(gcf, ['Unit ' num2str(i) '.jpg']);close all;
end
%exclude after assesement
exc_b=[36 45 48 119 122 170];

%exc_b=[];
%% 
%RED
cd('I:\neuropixels_artefact\onsets\red');
bs=1:200;
act=201:600;
fc=3;
time_onsetr=[];
for i=1:size(temp3,1)
trac=[];
trac=temp3(i,:);
[time_onsetr(i)] = onset_detection(trac,bs,act,fc);
%title(['unit ' num2str(i)]);saveas(gcf, ['Unit ' num2str(i) '.jpg']);close all;
end
%exclude after assesement
exc_r=[15 30 42 65 79 140 143 166 174 194 199 204 246 265];
%exc_r=[];
%% 
%ORANGE
cd('I:\neuropixels_artefact\onsets\orange');
bs=1:200;
act=201:600;
fc=3; 

time_onseto=[];
for i=1:size(temp2,1)
trac=[];
trac=temp2(i,:);
[time_onseto(i)] = onset_detection(trac,bs,act,fc);
%title(['unit ' num2str(i)]);saveas(gcf, ['Unit ' num2str(i) '.jpg']);close all;
end
%exclude after assesement
exc_o=[2 21 40 82 108 111 143 230];
%exc_o=[];
%% set manually excluded units to NAN
time_blue=time_onsetb;
time_orange=time_onseto;
time_red=time_onsetr;

time_blue(exc_b)=NaN;
time_orange(exc_o)=NaN;
time_red(exc_r)=NaN;


%% Violin plot time onsets at max 15 mW
g1=[];g2=[];g3=[];
p1=[];p2=[];p3=[];
p1=time_blue;p2=time_orange;p3=time_red;
p1(isnan(p1))=[];
p2(isnan(p2))=[];
p3(isnan(p3))=[];
color_id={[1 0 0],([255 165 0]/256),([0 191 255]/256)};
par=[];par=[p1 p2 p3]';
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
gro=[];gro=[g1 g2 g3]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 250]);set(gcf,'color','w');
violins = violinplot(par, gro,'ViolinColor',[([0 191 255]/256);([255 165 0]/256);[1 0 0]],'ShowMean', false,'ShowMedian', true,'ViolinAlpha',0,'MedianColor',[0 0 0],'BoxColor',[0 0 0]...
,'MarkerSize',10,'MedianMarkerSize',20);
box off;
ylim([0 400]);yticks([0:80:400])
xlim([0 4]);ylabel('Onset time (ms)');set(gca,'FontSize',12);set(gca,'TickDir','out');
h = gca;h.XAxis.Visible = 'off';
set(gca,'xtick',[]);
text(1.4,-60,'437 nm','Color',blue_c,'FontSize',11);
text(2.4,-60,'594 nm','Color',orange_c,'FontSize',11);
text(3.4,-60,'637 nm','Color',red_c,'FontSize',11);
line([0.9 3],[400 400],'Color','k');
line([1 2],[360 360],'Color','k');
t1=text(1.8,420,'***','FontSize',18,'Rotation',90);
t1=text(1.2,380,'***','FontSize',18,'Rotation',90)
title('15 mW','FontWeight','normal')
view([90 -90])
cd(save_folder);saveas(gcf, 'violin_latencies_15mW.pdf');
%display mean and SEM
disp('Time onset BLUE, ORANGE, RED')
[nanmean(p1) nanstd(p1)/sqrt(length(p1))]
[nanmean(p2) nanstd(p2)/sqrt(length(p2))]
[nanmean(p3) nanstd(p3)/sqrt(length(p3))]

%% Statisics for time onset
data_times=[];data_times=[p1' ;p2' ;p3'];
group_id=[];group_id=[g1'; g2'; g3'];
[p,tbl,stats]  = kruskalwallis(data_times,group_id);
figure;c = multcompare(stats);

%% LOWEST INTENSITY TIMING
%% Time onset BLUE, ORANGE, RED
p1=[];p2=[];p3=[];
temp1=[];temp2=[];temp3=[];
inten_i=1;
sess_i=1;

p1=zetas_blue(:,inten_i,sess_i);
p2=zetas_orange(:,inten_i,sess_i);
p3=zetas_red(:,inten_i,sess_i);

%only include <0.05 p values 
temp1=sdfs_blue(p1<0.05,:,inten_i,sess_i);
temp2=sdfs_orange(p2<0.05,:,inten_i,sess_i);
temp3=sdfs_red(p3<0.05,:,inten_i,sess_i);
%% ONSET TIMING lowest laser intensity
%BLUE
cd('I:\neuropixels_artefact\onsets\blue_low');
bs=1:200;
act=201:600;
fc=3;
time_onsetb=[];
for i=1:size(temp1,1)
trac=[];
trac=temp1(i,:);
[time_onsetb(i)] = onset_detection(trac,bs,act,fc);
%title(['unit ' num2str(i)]);saveas(gcf, ['Unit ' num2str(i) '.jpg']);close all;
end
exc_b=[];
%exclude after assesement
exc_b=[13 29 35 119 119 120];

%% 
%RED
cd('I:\neuropixels_artefact\onsets\red_low');
bs=1:200;
act=201:600;
fc=3;
time_onsetr=[];
for i=1:size(temp3,1)
trac=[];
trac=temp3(i,:);
[time_onsetr(i)] = onset_detection(trac,bs,act,fc);
%title(['unit ' num2str(i)]);saveas(gcf, ['Unit ' num2str(i) '.jpg']);close all;
end
%exclude after assesement
exc_r=[];
exc_r=[58 83 84 119 120 142 157 166 186 193 221 237 251];

%% 
%ORANGE
cd('I:\neuropixels_artefact\onsets\orange_low');
bs=1:200;
act=201:600;
fc=3; 

time_onseto=[];
for i=1:size(temp2,1)
trac=[];
trac=temp2(i,:);
[time_onseto(i)] = onset_detection(trac,bs,act,fc);
%title(['unit ' num2str(i)]);saveas(gcf, ['Unit ' num2str(i) '.jpg']);close all;
end
%exclude after assesement
exc_o=[];
exc_o=[7 40];

%% set manually excluded units to NAN
time_blue=time_onsetb;
time_orange=time_onseto;
time_red=time_onsetr;
time_blue(exc_b)=NaN;
time_orange(exc_o)=NaN;
time_red(exc_r)=NaN;
%% Violin plot time onsets at max 1 mW
g1=[];g2=[];g3=[];
p1=[];p2=[];p3=[];
p1=time_blue;p2=time_orange;p3=time_red;
p1(isnan(p1))=[];
p2(isnan(p2))=[];
p3(isnan(p3))=[];
color_id={[1 0 0],([255 165 0]/256),([0 191 255]/256)};
par=[];par=[p1 p2 p3]';
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
gro=[];gro=[g1 g2 g3]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 250]);set(gcf,'color','w');
violins = violinplot(par, gro,'ViolinColor',[([0 191 255]/256);([255 165 0]/256);[1 0 0]],'ShowMean', false,'ShowMedian', true,'ViolinAlpha',0,'MedianColor',[0 0 0],'BoxColor',[0 0 0]...
,'MarkerSize',10,'MedianMarkerSize',20);
box off;
ylim([0 400]);yticks([0:80:400])
xlim([0 4]);ylabel('Onset time (ms)');set(gca,'FontSize',12);set(gca,'TickDir','out');
h = gca;h.XAxis.Visible = 'off';
set(gca,'xtick',[]);
text(1.4,-60,'437 nm','Color',blue_c,'FontSize',11);
text(2.4,-60,'594 nm','Color',orange_c,'FontSize',11);
text(3.4,-60,'637 nm','Color',red_c,'FontSize',11);
line([1 3],[400 400],'Color','k');
line([2 3],[380 380],'Color','k');
text(2.15,400,'***','FontSize',18,'Rotation',90);
text(1.9,420,'***','FontSize',18,'Rotation',90)
title('1 mW','FontWeight','normal')
view([90 -90])
cd(save_folder);saveas(gcf, 'violin_latencies_1mW.pdf');
% display mean and SEM
disp('Time onset BLUE, ORANGE, RED')
[nanmean(p1) nanstd(p1)/sqrt(length(p1))]
[nanmean(p2) nanstd(p2)/sqrt(length(p2))]
[nanmean(p3) nanstd(p3)/sqrt(length(p3))]
%% Statisics for time onset
data_times=[];data_times=[p1' ;p2' ;p3'];
group_id=[];group_id=[g1'; g2'; g3'];
[p,tbl,stats]  = kruskalwallis(data_times,group_id);
figure;c = multcompare(stats);
%% AMBIENT LIGHT 
%% Plot delta firing rate dark vs light only the ones that were there in dark Fig2C
% blue 
a=[];a=find(zetas_blue(:,5,1)<0.05);
data=[];data=[MIs_blue(a,5,1) MIs_blue(a,5,2)];
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 900, 300]);set(gcf,'color','w');t=tiledlayout("horizontal",'TileSpacing','Compact');
 nexttile;
%xlim([-10 40]);ylim([-10 40]);xticks([-10:10:40]);yticks([-10:10:40])
xlim([-30 60]);ylim([-30 60]);xticks([-30:30:60]);yticks([-30:30:60])
 xlabel({'\Delta Firing rate (Hz)'; 'LIGHT (20 lux)'});ylabel({'DARK';'Delta Firing rate (Hz)'});
%text(50,500,['p=' num2str(round(p_ct,2))]);
%hold on;plot([0 max([xlim ylim])], [0 max([xlim ylim])], '--k');
rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '-';
line([0 0],[-30 60],'Color','k','LineStyle','--');
line([-30 60],[0 0],'Color','k','LineStyle','--');
hold on;
 s1=scatter(MIs_blue(a,5,2),MIs_blue(a,5,1),40,blue_c,'filled','o');axis square;hold on;
 s1.MarkerEdgeColor = [0 0 0];axis square;
set(gca,'FontSize',11);set(gca,'TickDir','out'); %title('Excitation','Color','r','FontWeight','normal');
 offsetAxes;xtickangle(0);text(30,60,'437 nm','Color',blue_c,'FontSize',11);
  %stats
nanmean(data)
nanstd(data)/sqrt(length(data))
length(data)
[h,p] = adtest(data(:,1))
[h,p]=signrank(data(:,1),data(:,2))
%orange
nexttile
a=[];a=find(zetas_orange(:,5,1)<0.05);
data=[];data=[MIs_orange(a,5,1) MIs_orange(a,5,2)];
 %xlim([-30 40]);ylim([-30 40]);xticks([-30:10:40]);yticks([-30:10:40])
 xlim([-30 60]);ylim([-30 60]);xticks([-30:30:60]);yticks([-30:30:60])
 rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '-';
line([0 0],[-30 60],'Color','k','LineStyle','--');
line([-30 60],[0 0],'Color','k','LineStyle','--');
hold on;
 s1=scatter(MIs_orange(a,5,2),MIs_orange(a,5,1),40,orange_c,'filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
 s1.MarkerEdgeColor = [0 0 0];axis square;
;%xticks([0:250:500]);yticks([0:250:500])
 xlabel({'\Delta Firing rate (Hz)'; 'LIGHT (20 lux)'});ylabel({'DARK';'Delta Firing rate (Hz)'});
 set(gca,'FontSize',11);set(gca,'TickDir','out');
title('15 mW','FontWeight','normal')
offsetAxes;xtickangle(0);h = gca;h.YAxis.Visible = 'off';text(30,60,'594 nm','Color',orange_c,'FontSize',11);
 %stats
nanmean(data)
nanstd(data)/sqrt(length(data))
%[r p] = signrank(data(:,1),data(:,2))
length(data)
[h,p] = adtest(data(:,1))
[h,p]=signrank(data(:,1),data(:,2))
 %red
 nexttile
a=[];a=find(zetas_red(:,5,1)<0.05);
data=[];data=[MIs_red(a,5,1) MIs_red(a,5,2)];

xlim([-30 60]);ylim([-30 60]);xticks([-30:30:60]);yticks([-30:30:60])
 rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '-';
line([0 0],[-30 60],'Color','k','LineStyle','--');
line([-30 60],[0 0],'Color','k','LineStyle','--');
hold on;
 s1=scatter(MIs_red(a,5,2),MIs_red(a,5,1),40,red_c,'filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
 s1.MarkerEdgeColor = [0 0 0];axis square;
 xlabel({'\Delta Firing rate (Hz)'; 'LIGHT (20 lux)'});ylabel({'DARK';'Delta Firing rate (Hz)'});
 set(gca,'FontSize',11);set(gca,'TickDir','out');
 offsetAxes;h = gca;h.YAxis.Visible = 'off';text(30,60,'637 nm','Color',red_c,'FontSize',11);
 %stats
nanmean(data)
nanstd(data)/sqrt(length(data))
length(data)
[h,p] = adtest(data(:,1))
[h,p]=signrank(data(:,1),data(:,2))
% save
cd(save_folder);saveas(gcf, 'comparison units from dark to light at 15 mW.pdf');
 
%% Baseline vs light response PERMUTATION test 20 lux test highest intensity
%change session and intesity accordingly 
sess_i=2;intens_i=5;
baser=1:200;
respr=201:700;
temp_sdf={sdfs_blue,sdfs_orange,sdfs_red};
tmp_zeta={zetas_blue,zetas_orange,zetas_red};
p_all=[];p=[];
for k=1:10
for i=1:3
    sdf=[];sdf=temp_sdf{i};
    zta=[];zta=tmp_zeta{i};trz=[];tr1=[];tr2=[];
    tr1=nanmean(sdf(:,baser,intens_i,sess_i),2);
    tr2=nanmean(sdf(:,respr,intens_i,sess_i),2);
    trz=zta(:,intens_i,sess_i);
 %[r(i) p(i)]=signrank(tr1(~isnan(trz)), tr2(~isnan(trz)));
 %[p(i), observeddifference, effectsize] = permutationTest(tr1(~isnan(trz)), tr2(~isnan(trz)), 1000);
 [p(i), observeddifference, effectsize] = permutationTest(tr1, tr2, 1000);

mb(i)=nanmean(tr1(~isnan(trz)));
mr(i)=nanmean(tr2(~isnan(trz)));
mb_sem(i)=nanstd(tr1(~isnan(trz)))/sqrt(length(tr1(~isnan(trz))));
mr_sem(i)=nanstd(tr2(~isnan(trz)))/sqrt(length(tr2(~isnan(trz))));

    %[h(i) p] = adtest(nanmean(sdf(:,respr,intens_i,sess_i),2));
    %[r(i) p(i)]=signrank(nanmean(sdf(:,baser,intens_i,sess_i),2), nanmean(sdf(:,respr,intens_i,sess_i),2));
  % [r(i) p(i)]=ranksum(allbase(:), nanmean(sdf(:,respr,intens_i,sess_i),2));
end
p_all(k,:)=p;
end
mb
mb_sem
mr
mr_sem
%% Baseline vs light response PERMUTATION test 20 lux test second lowest intensity
%change session and intesity accordingly 
sess_i=2;intens_i=2;
baser=1:200;
respr=201:700;
temp_sdf={sdfs_blue,sdfs_orange,sdfs_red};
tmp_zeta={zetas_blue,zetas_orange,zetas_red};
p_all=[];p=[];
for k=1:10
for i=1:3
    sdf=[];sdf=temp_sdf{i};
    zta=[];zta=tmp_zeta{i};trz=[];tr1=[];tr2=[];
    tr1=nanmean(sdf(:,baser,intens_i,sess_i),2);
    tr2=nanmean(sdf(:,respr,intens_i,sess_i),2);
    trz=zta(:,intens_i,sess_i);
 %[r(i) p(i)]=signrank(tr1(~isnan(trz)), tr2(~isnan(trz)));
 %[p(i), observeddifference, effectsize] = permutationTest(tr1(~isnan(trz)), tr2(~isnan(trz)), 1000);
 [p(i), observeddifference, effectsize] = permutationTest(tr1, tr2, 1000);

mb(i)=nanmean(tr1(~isnan(trz)));
mr(i)=nanmean(tr2(~isnan(trz)));
mb_sem(i)=nanstd(tr1(~isnan(trz)))/sqrt(length(tr1(~isnan(trz))));
mr_sem(i)=nanstd(tr2(~isnan(trz)))/sqrt(length(tr2(~isnan(trz))));

    %[h(i) p] = adtest(nanmean(sdf(:,respr,intens_i,sess_i),2));
    %[r(i) p(i)]=signrank(nanmean(sdf(:,baser,intens_i,sess_i),2), nanmean(sdf(:,respr,intens_i,sess_i),2));
  % [r(i) p(i)]=ranksum(allbase(:), nanmean(sdf(:,respr,intens_i,sess_i),2));
end
p_all(k,:)=p;
end
mb
mb_sem
mr
mr_sem
%% Table with all percentages 
    adder=1;
data_b=[];data_o=[];data_r=[];data_ba=[];data_oa=[];data_ra=[];data_bn=[];data_on=[];data_rn=[];data_bna=[];data_ona=[];data_rna=[];
zeta_temp1=[];zeta_temp2=[];zeta_temp3=[];zeta_temp4=[];zeta_temp5=[];zeta_temp6=[];
zeta_temp1=zetas_blue;zeta_temp2=zetas_orange;zeta_temp3=zetas_red;
zeta_temp4=zetasn_blue;zeta_temp5=zetasn_orange;zeta_temp6=zetasn_red;

for sess_i=1:4;

temp1=[];temp2=[];temp3=[];temp4=[];temp5=[];temp6=[];
temp1=(sum(zeta_temp1(:,:,sess_i)<0.05)/length(zeta_temp1))*100;
temp2=(sum(zeta_temp2(:,:,sess_i)<0.05)/length(zeta_temp2))*100;
temp3=(sum(zeta_temp3(:,:,sess_i)<0.05)/length(zeta_temp3))*100;
temp4=(sum(zeta_temp4(:,:,sess_i)<0.05)/length(zeta_temp4))*100;
temp5=(sum(zeta_temp5(:,:,sess_i)<0.05)/length(zeta_temp5))*100;
temp6=(sum(zeta_temp6(:,:,sess_i)<0.05)/length(zeta_temp6))*100;

data_b(adder,:)=[round(temp1,2)];
data_o(adder,:)=[round(temp2,2)];
data_r(adder,:)=[round(temp3,2)];
data_bn(adder,:)=[round(temp4,2)];
data_on(adder,:)=[round(temp5,2)];
data_rn(adder,:)=[round(temp6,2)];
data_ba(adder,:)=[sum(zeta_temp1(:,:,sess_i)<0.05)];
data_oa(adder,:)=[sum(zeta_temp2(:,:,sess_i)<0.05)];
data_ra(adder,:)=[sum(zeta_temp3(:,:,sess_i)<0.05)];
data_bna(adder,:)=[sum(zeta_temp4(:,:,sess_i)<0.05)];
data_ona(adder,:)=[sum(zeta_temp5(:,:,sess_i)<0.05)];
data_rna(adder,:)=[sum(zeta_temp6(:,:,sess_i)<0.05)];
adder=adder+1;

end
%% Plot heatmaps with percentage and control ZETA for Fig2D
d_all={[data_b ;nanmean(data_bn)],[data_o ;nanmean(data_on)],[data_r ;nanmean(data_rn)]};
color_v=[blue_c ;orange_c ;red_c];
for i=1:3
    data_m=[];
data_m=d_all{i};
data_m=[data_m(1,:) ;[40 40 40 40 40] ;data_m(2:end,:)]

fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400+i*300, 200, 300, 200]);set(gcf,'color','w');
 imagesc(data_m);
cmap(color_v(i,:),100,20,0.5);h=colorbar;caxis([0 40]);
for k=1:5  
text(k-0.25,1,num2str(round(data_m(1,k),1)));
hold on
end
for k=1:5
text(k-0.25,3,num2str(round(data_m(3,k),1)),'Color','w');
end
for k=1:5
text(k-0.25,4,num2str(round(data_m(4,k),1)),'Color','w');
end
for k=1:5
text(k-0.25,5,num2str(round(data_m(5,k),1)),'Color','w');
end
for k=1:5
text(k-0.25,6,num2str(round(data_m(6,k),1)),'Color','w');
end
h.Label.String = {'(%)'};
ylabel('Ambient light (lux)');xlabel('Laser power (mW)');
xticklabels({'1','2.5','5','10','15'});yticks([1:1:6]);yticklabels({'dark','','20','40','80','control'});set(gca,'FontSize',11);box off;
line([0.5 5.5],[5.5 5.5],'Color','w','LineStyle','--');
cd(save_folder);saveas(gcf, [num2str(i) ' table_light.pdf']);
end
%% test against zeta noise
sess_i=4;
testa1=data_ra;
testa2=nanmean(data_rna);p=[];
for i=1:5
    x=[];x=[testa1(sess_i,i),length(MIs_blue)-testa1(sess_i,i);round(testa2(i)),length(MIs_blue)-round(testa2(i))];
    %[h,p(i),stats] = fishertest(x,'Tail','right');
    [h,p(i),stats] = fishertest(x); 
end
p

%% SUPPLEMENTARY MATERIAL
%% heatmap time display BLUE DARK Fig S1
sess_i=1;
color_v=blue_c;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 1000, 300]);set(gcf,'color','w');
t = tiledlayout(1,5,'TileSpacing','Compact');
display_s=1:700;
for i=1:5
temp1=[];
temp1=sdfs_blue(:,display_s,i,sess_i);
max(max(temp1))
nexttile;
imagesc(temp1);
cc=colormap(cmap(color_v,100,5,2));caxis([0 100]);
line([200 200],[0 max(get(gca,'ylim'))],'Color','w','LineStyle','--');
h=gca;h.YAxis.Visible = 'off'; 
box off;set(gca,'TickDir','out');
yticklabels({''})
xticks([1:200:700])
xticklabels({'-0.2','0','0.2','0.4','0.6','0.4','0.5'})
xlabel('Time (s)');set(gca,'FontSize',11);
title([num2str(l_power(i)) ' mW'],'FontWeight','normal');
if i==1
    t=text(-50,length(temp1)/1.5,'Units','FontSize',12);
     set(t,'Rotation',90);
end
end
h=colorbar;
h.Label.String = 'Firing rate (Hz)';set(gca,'FontSize',11);
cd(save_folder);saveas(gcf, 'units_blue.pdf');
%% heatmap time display ORANGE DARK Fig S1 2nd row
sess_i=1;
color_v=orange_c;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 1000, 300]);set(gcf,'color','w');
t = tiledlayout(1,5,'TileSpacing','Compact');
display_s=1:700;
for i=1:5
temp1=[];
temp1=sdfs_orange(:,display_s,i,sess_i);
max(max(temp1))
nexttile;
imagesc(temp1);
cc=colormap(cmap(color_v,100,5,2));caxis([0 100]);
line([200 200],[0 max(get(gca,'ylim'))],'Color','w','LineStyle','--');
h=gca;h.YAxis.Visible = 'off'; 
box off;set(gca,'TickDir','out');
yticklabels({''})
xticks([1:200:700])
xticklabels({'-0.2','0','0.2','0.4','0.6','0.4','0.5'})
xlabel('Time (s)');set(gca,'FontSize',11);
title([num2str(l_power(i)) ' mW'],'FontWeight','normal');
if i==1
    t=text(-50,length(temp1)/1.5,'Units','FontSize',12);
     set(t,'Rotation',90);
end
end
h=colorbar;
h.Label.String = 'Firing rate (Hz)';set(gca,'FontSize',11);
cd(save_folder);saveas(gcf, 'units_orange.pdf');
%% heatmap time display RED DARK 3rd row
sess_i=1;
color_v=red_c;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 1000, 300]);set(gcf,'color','w');
t = tiledlayout(1,5,'TileSpacing','Compact');
display_s=1:700;
for i=1:5
temp1=[];
temp1=sdfs_red(:,display_s,i,sess_i);
max(max(temp1))
nexttile;
imagesc(temp1);
cc=colormap(cmap(color_v,100,5,2));caxis([0 100]);
line([200 200],[0 max(get(gca,'ylim'))],'Color','w','LineStyle','--');
h=gca;h.YAxis.Visible = 'off'; 
box off;set(gca,'TickDir','out');
yticklabels({''})
xticks([1:200:700])
xticklabels({'-0.2','0','0.2','0.4','0.6','0.4','0.5'})
xlabel('Time (s)');set(gca,'FontSize',11);
title([num2str(l_power(i)) ' mW'],'FontWeight','normal');
if i==1
    t=text(-50,length(temp1)/1.5,'Units','FontSize',12);
     set(t,'Rotation',90);
end
end
h=colorbar;
h.Label.String = 'Firing rate (Hz)';set(gca,'FontSize',11);
cd(save_folder);saveas(gcf, 'units_red.pdf');

%% heatmap time display BLUE 20 LUX FIG S2A
sess_i=2;
color_v=blue_c;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 1000, 300]);set(gcf,'color','w');
t = tiledlayout(1,5,'TileSpacing','Compact');
display_s=1:700;
for i=1:5
temp1=[];
temp1=sdfs_blue(:,display_s,i,sess_i);
max(max(temp1))
nexttile;
imagesc(temp1);
cc=colormap(cmap(color_v,100,5,2));caxis([0 100]);
line([200 200],[0 max(get(gca,'ylim'))],'Color','w','LineStyle','--');
h=gca;h.YAxis.Visible = 'off'; 
box off;set(gca,'TickDir','out');
yticklabels({''})
xticks([1:200:700])
xticklabels({'-0.2','0','0.2','0.4','0.6','0.4','0.5'})
xlabel('Time (s)');set(gca,'FontSize',11);
title([num2str(l_power(i)) ' mW'],'FontWeight','normal');
if i==1
    t=text(-50,length(temp1)/1.5,'Units','FontSize',12);
     set(t,'Rotation',90);
end
end
h=colorbar;
h.Label.String = 'Firing rate (Hz)';set(gca,'FontSize',11);
cd(save_folder);saveas(gcf, 'units_blue_20lux.pdf');
%% heatmap time display ORANGE 20 lux FIG S2A 2nd row
sess_i=2;
color_v=orange_c;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 1000, 300]);set(gcf,'color','w');
t = tiledlayout(1,5,'TileSpacing','Compact');
display_s=1:700;
for i=1:5
temp1=[];
temp1=sdfs_orange(:,display_s,i,sess_i);
max(max(temp1))
nexttile;
imagesc(temp1);
cc=colormap(cmap(color_v,100,5,2));caxis([0 100]);
line([200 200],[0 max(get(gca,'ylim'))],'Color','w','LineStyle','--');
h=gca;h.YAxis.Visible = 'off'; 
box off;set(gca,'TickDir','out');
yticklabels({''})
xticks([1:200:700])
xticklabels({'-0.2','0','0.2','0.4','0.6','0.4','0.5'})
xlabel('Time (s)');set(gca,'FontSize',11);
title([num2str(l_power(i)) ' mW'],'FontWeight','normal');
if i==1
    t=text(-50,length(temp1)/1.5,'Units','FontSize',12);
     set(t,'Rotation',90);
end
end
h=colorbar;
h.Label.String = 'Firing rate (Hz)';set(gca,'FontSize',11);
cd(save_folder);saveas(gcf, 'units_orange_20lux.pdf');
%% heatmap time display RED 20 lux FIG S2A 3rd row
sess_i=2;
color_v=red_c;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 1000, 300]);set(gcf,'color','w');
t = tiledlayout(1,5,'TileSpacing','Compact');
display_s=1:700;
for i=1:5
temp1=[];
temp1=sdfs_red(:,display_s,i,sess_i);
max(max(temp1))
nexttile;
imagesc(temp1);
cc=colormap(cmap(color_v,100,5,2));caxis([0 100]);
line([200 200],[0 max(get(gca,'ylim'))],'Color','w','LineStyle','--');
h=gca;h.YAxis.Visible = 'off'; 
box off;set(gca,'TickDir','out');
yticklabels({''})
xticks([1:200:700])
xticklabels({'-0.2','0','0.2','0.4','0.6','0.4','0.5'})
xlabel('Time (s)');set(gca,'FontSize',11);
title([num2str(l_power(i)) ' mW'],'FontWeight','normal');
if i==1
    t=text(-50,length(temp1)/1.5,'Units','FontSize',12);
     set(t,'Rotation',90);
end
end
h=colorbar;
h.Label.String = 'Firing rate (Hz)';set(gca,'FontSize',11);
cd(save_folder);saveas(gcf, 'units_red_20lux.pdf');



%% Percentage modified plot for 20 lux and darknes in one FIG S2B

%single or multi
zeta_temp1=[];zeta_temp2=[];zeta_temp3=[];temp1a=[];temp2a=[];temp3a=[];
zeta_temp1=zetas_blue;zeta_temp2=zetas_orange;zeta_temp3=zetas_red;
temp1=[];temp2=[];temp3=[];temp4=[];temp5=[];temp6=[];

temp1=(sum(zeta_temp1(:,:,1)<0.05)/length(zeta_temp1))*100;
temp2=(sum(zeta_temp2(:,:,1)<0.05)/length(zeta_temp2))*100;
temp3=(sum(zeta_temp3(:,:,1)<0.05)/length(zeta_temp3))*100;

temp4=(sum(zeta_temp1(:,:,2)<0.05)/length(zeta_temp1))*100;
temp5=(sum(zeta_temp2(:,:,2)<0.05)/length(zeta_temp2))*100;
temp6=(sum(zeta_temp3(:,:,2)<0.05)/length(zeta_temp3))*100;
%PLOT
fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [900, 400 ,800, 300]);
b1=bar([temp1(1) temp2(1) temp3(1) 0 temp4(1) temp5(1) temp6(1) ;temp1(2) temp2(2) temp3(2) 0 temp4(2) temp5(2) temp6(2) ;temp1(3) temp2(3) temp3(3) 0 temp4(3) temp5(3) temp6(3) ;...
    temp1(4) temp2(4) temp3(4) 0 temp4(4) temp5(4) temp6(4) ;temp1(5) temp2(5) temp3(5) 0 temp4(5) temp5(5) temp6(5)]);
set(b1,'ShowBaseLine','off');b1(1).FaceColor=blue_c;b1(2).FaceColor=orange_c;b1(3).FaceColor=red_c;b1(4).FaceColor=[1 1 1];b1(5).FaceColor=blue_c;b1(6).FaceColor=orange_c;b1(7).FaceColor=red_c;
b1(1).LineWidth=1;b1(2).LineWidth=1;b1(3).LineWidth=1;b1(4).LineWidth=0.0001;b1(5).LineWidth=1;b1(6).LineWidth=1;b1(7).LineWidth=1;
%b1(1).EdgeColor=[0 0 0];b2(1).EdgeColor=[0 0 0];b3(1).EdgeColor=[0 0 0];b4(1).EdgeColor=[0.5 0.5 0.5];b5(1).EdgeColor=[0.5 0.5 0.5];b6(1).EdgeColor=[0.5 0.5 0.5];b7(1).EdgeColor=[0.5 0.5 0.5];
ylim([0 50]);

box off;set(gca,'TickDir','out');ylabel('Percentage modulated units (%)');xlabel('Fiber laser power (mW)');
offsetAxes;hold on;text(6.5,max(get(gca,'ylim')),'637 nm','Color','r','FontWeight','bold','FontSize',12);
hold on;text(6.5,max(get(gca,'ylim')-4),'594 nm','Color',([255 165 0]/256),'FontWeight','bold','FontSize',12);
hold on;text(6.5,max(get(gca,'ylim')-8),'473 nm ','Color',([0 191 255]/256),'FontWeight','bold','FontSize',12);xticklabels({'1','2.5','5','10','15'});set(gca,'FontSize',12);
text(0.55,45,'DARK','FontSize',12);text(1.05,45,'LIGHT','FontSize',11,'Color',[0.5 0.5 0.5]);
text(1.55,45,'DARK','FontSize',12);text(2.05,45,'LIGHT','FontSize',11,'Color',[0.5 0.5 0.5]);
text(2.55,45,'DARK','FontSize',12);text(3.05,45,'LIGHT','FontSize',11,'Color',[0.5 0.5 0.5]);
text(3.55,45,'DARK','FontSize',12);text(4.05,45,'LIGHT','FontSize',11,'Color',[0.5 0.5 0.5]);
text(4.55,45,'DARK','FontSize',12);text(5.05,45,'LIGHT','FontSize',11,'Color',[0.5 0.5 0.5]);
%save
cd(save_folder);saveas(gcf, 'percentage_light_compare_dark_light.pdf');



%%%%%%%%%%%%%%REVIEWS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Referee question about latencies of cortical layers 
for i=1:4
 %SUA
l23_b{i}=find(blue_depth==821 & blue_a==i);
l4_b{i}=find(blue_depth==721 & blue_a==i);
l5_b{i}=find(blue_depth==778 & blue_a==i);
l6_b{i}=find(blue_depth==33 & blue_a==i | blue_depth==305 & blue_a==i);

l23_o{i}=find(orange_depth==821 & orange_a==i);
l4_o{i}=find(orange_depth==721 & orange_a==i);
l5_o{i}=find(orange_depth==778 & orange_a==i);
l6_o{i}=find(orange_depth==33 & orange_a==i | orange_depth==305 & orange_a==i);

l23_r{i}=find(red_depth==821 & red_a==i);
l4_r{i}=find(red_depth==721 & red_a==i);
l5_r{i}=find(red_depth==778 & red_a==i);
l6_r{i}=find(red_depth==33 & red_a==i | red_depth==305 & red_a==i);

%MUA
l23_bm{i}=find(blue_depthm==821 & blue_am==i);
l4_bm{i}=find(blue_depthm==721 & blue_am==i);
l5_bm{i}=find(blue_depthm==778 & blue_am==i);
l6_bm{i}=find(blue_depthm==33 & blue_am==i | blue_depthm==305 & blue_am==i);

l23_om{i}=find(orange_depthm==821 & orange_am==i);
l4_om{i}=find(orange_depthm==721 & orange_am==i);
l5_om{i}=find(orange_depthm==778 & orange_am==i);
l6_om{i}=find(orange_depthm==33 & orange_am==i | orange_depthm==305 & orange_am==i);

l23_rm{i}=find(red_depthm==821 & red_am==i);
l4_rm{i}=find(red_depthm==721 & red_am==i);
l5_rm{i}=find(red_depthm==778 & red_am==i);
l6_rm{i}=find(red_depthm==33 & red_am==i | red_depthm==305 & red_am==i);
end
%% Plot average per animal ANIMAL 4 blue 
sdf_tm=[];
sdf_tm=temp1(:,1:700)-nanmean(temp1(:,1:199),2);
ani=4;
color1='m';color2='m';color3='g';color4='r';
stim_on=20;

tx1=[];tx2=[];tx3=[];tx4=[];
tx1=sdf_tm(l5_b{ani},:);
tx2=sdf_tm(l5_b{ani},:);
tx3=sdf_tm(l6_b{ani},:);
tx4=sdf_tm(l4_b{ani},:);
average_errorbar_sdf_layer(tx1,tx2,tx3,tx4,stim_on,color1,color2,color3,color4)
title('Blue illumination','FontWeight','normal')
text(320,28,'L4','Color',color4);
text(320,24,'L5','Color',color1);
text(320,20,'L6','Color',color3);
ylim([-5 30]);xlim([150 350]);
hold on;line([stim_on*10, stim_on*10], [0 max(get(gca, 'ylim'))-0.5], 'color', 'b', 'linestyle', '--')
%offsetAxes;
%% Plot average per animal ANIMAL 4 orange
sdf_tm=[];
sdf_tm=temp2(:,1:700)-nanmean(temp2(:,1:199),2);
ani=4;
color1='m';color2='m';color3='g';color4='r';
stim_on=20;

tx1=[];tx2=[];tx3=[];tx4=[];
tx1=sdf_tm(l5_o{ani},:);
tx2=sdf_tm(l5_o{ani},:);
tx3=sdf_tm(l6_o{ani},:);
tx4=sdf_tm(l4_o{ani},:);
average_errorbar_sdf_layer(tx1,tx2,tx3,tx4,stim_on,color1,color2,color3,color4)
title('Orange illumination animal 4','FontWeight','normal')
text(600,28,'L4','Color',color4);
text(600,24,'L5','Color',color1);
text(600,20,'L6','Color',color3);
ylim([-2 33]);xlim([200 300])
offsetAxes;
%% Plot average per animal ANIMAL 4 RED
sdf_tm=[];
sdf_tm=temp3(:,1:700)-nanmean(temp3(:,1:199),2);
ani=4;
color1='m';color2='m';color3='g';color4='r';
stim_on=20;

tx1=[];tx2=[];tx3=[];tx4=[];
tx1=sdf_tm(l5_r{ani},:);
tx2=sdf_tm(l5_r{ani},:);
tx3=sdf_tm(l6_r{ani},:);
tx4=sdf_tm(l4_r{ani},:);
average_errorbar_sdf_layer(tx1,tx2,tx3,tx4,stim_on,color1,color2,color3,color4)
title('Red illumination animal 4','FontWeight','normal')
text(600,28,'L4','Color',color4);
text(600,24,'L5','Color',color1);
text(600,20,'L6','Color',color3);
ylim([-2 30]);xlim([200 300])
offsetAxes;





%% ATTACHED FUNCTIONS

function [sdf_sua sdf_mua para_sua para_mua zeta_sua zeta_mua zeta_noise_sua depth_id depth_nr depth_idm depth_nrm] = across_animals(neuropixel_data,sess_i,l_inten)

depth_idm=[];depth_nrm=[];depth_nr=[];depth_id=[];sdf_sua=[];sdf_mua=[];para_sua=[];para_mua=[];zeta_sua=[];zeta_mua=[];zeta_noise_sua=[];
for i=1:length(neuropixel_data)
sdf_sua=[sdf_sua ;neuropixel_data(i).session_aligned(sess_i).sdf_cortex(:,:,l_inten)];
sdf_mua=[sdf_mua ;neuropixel_data(i).session_aligned(sess_i).sdf_mua(:,:,l_inten)];
para_sua=[para_sua;neuropixel_data(i).session_aligned(sess_i).param_cortex(:,:,l_inten)];
para_mua=[para_mua;neuropixel_data(i).session_aligned(sess_i).param_mua(:,:,l_inten)];
zeta_sua=[zeta_sua;neuropixel_data(i).session_aligned(sess_i).zeta_cortex(:,l_inten)];
zeta_noise_sua=[zeta_noise_sua;neuropixel_data(i).session_aligned(sess_i).zeta_cortex_noise(:,l_inten)];
zeta_mua=[zeta_mua;neuropixel_data(i).session_aligned(sess_i).zeta_mua(:,l_inten)];
depth_id=[depth_id ;neuropixel_data(i).cortex_id];
depth_nr=[depth_nr; cellfun(@(v)v(1),neuropixel_data(i).cortex_depth)];
depth_idm=[depth_idm ;neuropixel_data(i).cortexmua_id];
depth_nrm=[depth_nrm; cellfun(@(v)v(1),neuropixel_data(i).cortexmua_depth)];
end
end


function [time_onset] = onset_detection(trac,bs,act,fc)
filter_used='gaussian';
filter_trac=smoothdata(trac,filter_used,50);
diff_trac=diff(filter_trac);
bs_diff=diff_trac(bs);
act_diff=diff_trac(act);
st_bs=fc*std(bs_diff);
%trace_smooth=trac(act);

if st_bs<0.1
 
    [st_bs, locs] = findpeaks(act_diff);
end
try
temp=find(act_diff>=st_bs(1));
catch
    temp=NaN;
end


 % fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [900, 400 ,400, 300]);
 % plot(trac,'Color','k','LineWidth',1.5);box off;hold on;
 % line([bs(end) bs(end)],[0 max(get(gca,'ylim'))],'Color','b');
 % offsetAxes;set(gca,'TickDir','out');


try
time_onset=temp(1);

line([time_onset+bs(end) time_onset+bs(end)],[0 max(get(gca,'ylim'))],'Color','r');
text(time_onset+bs(end),max(get(gca,'ylim')),[num2str(time_onset) ' ms']);
catch
    time_onset=NaN;
end



end

%% 
function average_errorbar_sdf(temp1,temp2,temp3,sess_i,inten_i,scalef,l_power,stim_on,blue_c,orange_c,red_c,lux_m)
%PLOT
 fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [900, 400 ,350, 225]);
 s2=shadedErrorBar(1:size(temp1,2),nanmean(temp1,1),nanstd(temp1,1)/sqrt(size(temp1,1)),'lineProps','k');hold on
 s1=shadedErrorBar(1:size(temp2,2),nanmean(temp2,1),nanstd(temp2,1)/sqrt(size(temp2,1)),'lineProps','k');hold on
 s3=shadedErrorBar(1:size(temp3,2),nanmean(temp3,1),nanstd(temp3,1)/sqrt(size(temp3,1)),'lineProps','k');hold on
set(gca,'XTickLabel',[]);hold on;box off
xticklabels({'-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'});ylabel('Firing rate (Hz)');xlabel('Time (s)');xticks([0 100 200 300 400 500 600 700]);
ylim([0 10]);%xlim([20 500]);
hold on;line([stim_on*10, stim_on*10], [0 max(get(gca, 'ylim'))-0.5], 'color', 'k', 'linestyle', '--');xlim([100 size(temp1,2)-500]);set(gca,'FontSize',12);
s3.mainLine.Color=red_c;s1.mainLine.LineWidth=2;
s1.mainLine.Color=orange_c;s2.mainLine.LineWidth=2;
s2.mainLine.Color=blue_c;s3.mainLine.LineWidth=2;
%plot([400 420],[max(get(gca,'ylim')) max(get(gca,'ylim'))],'Color','r');
hold on;text(423,max(get(gca,'ylim')),'637 nm','Color','r','FontWeight','bold','FontSize',12);
%plot([400 420],[max(get(gca,'ylim'))-0.5*scalef max(get(gca,'ylim'))-0.5*scalef],'Color',([255 165 0]/256));
 hold on;text(423,max(get(gca,'ylim')-0.75*scalef),'594 nm','Color',([255 165 0]/256),'FontWeight','bold','FontSize',12);
%plot([400 420],[max(get(gca,'ylim'))-1*scalef max(get(gca,'ylim'))-1*scalef],'Color',([0 191 255]/256));
hold on;text(423,max(get(gca,'ylim')-1.5*scalef),'473 nm ','Color',([0 191 255]/256),'FontWeight','bold','FontSize',12);
hold on;text(200,max(get(gca,'ylim')),'laser ON ','Color','k','FontWeight','normal','FontSize',12);
set(gca,'FontSize',12);title([' ' num2str(l_power(inten_i)) ' mW ' lux_m ' lux'],'FontWeight','normal');set(gca,'TickDir','out');
offsetAxes;h = gca;h.XAxis.Visible = 'off';  line([450 500],[0 0],'Color','k');text(450,-0.25,'50 ms');
end

function average_errorbar_sdf_layer(temp1,temp2,temp3,temp4,stim_on,color1,color2,color3,color4)
%PLOT
 fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [900, 400 ,500, 225]);
 s1=shadedErrorBar(1:size(temp1,2),nanmean(temp1,1),nanstd(temp1,1)/sqrt(size(temp1,1)),'lineProps','k','patchSaturation',0.05);hold on
 s2=shadedErrorBar(1:size(temp2,2),nanmean(temp2,1),nanstd(temp2,1)/sqrt(size(temp2,1)),'lineProps','k','patchSaturation',0.05);hold on
 s3=shadedErrorBar(1:size(temp3,2),nanmean(temp3,1),nanstd(temp3,1)/sqrt(size(temp3,1)),'lineProps','k','patchSaturation',0.05);hold on
 s4=shadedErrorBar(1:size(temp4,2),nanmean(temp4,1),nanstd(temp4,1)/sqrt(size(temp4,1)),'lineProps','k','patchSaturation',0.05);hold on
set(gca,'XTickLabel',[]);hold on;box off
xticklabels({'-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'});ylabel('Firing rate (Hz)');xlabel('Time (s)');xticks([0 100 200 300 400 500 600 700]);
%ylim([0 10]);%xlim([20 500]);
hold on;line([stim_on*10, stim_on*10], [0 max(get(gca, 'ylim'))-0.5], 'color', 'k', 'linestyle', '--');%xlim([100 size(temp1,2)-500]);
set(gca,'FontSize',12);
s1.mainLine.Color=color1;s1.mainLine.LineWidth=2;
s2.mainLine.Color=color2;s2.mainLine.LineWidth=2;
s3.mainLine.Color=color3;s3.mainLine.LineWidth=2;
s4.mainLine.Color=color4;s4.mainLine.LineWidth=2;
set(gca,'TickDir','out');

end


%backup 

% %% readout peak 
% sdf_tm=[];
% sdf_tm=temp1(:,1:700)-nanmean(temp1(:,1:199),2);
% ani=4;
% color1='m';color2='m';color3='g';color4='r';
% stim_on=20;
% 
% tx1=[];tx2=[];tx3=[];tx4=[];
% tx1=sdf_tm(l5_b{ani},:);
% tx2=sdf_tm(l5_b{ani},:);
% tx3=sdf_tm(l6_b{ani},:);
% tx4=sdf_tm(l4_b{ani},:);
% %% 
% 
% for i=1:size(tx4,1)
% figure;plot(tx4(i,:))
% [x(i) y(i)] = ginput(1);
% close all;
% end
% %% L5
% for i=1:size(tx1,1)
% figure;plot(tx1(i,:))
% [x5(i) y5(i)] = ginput(1);
% close all;
% end
% %% L6
% for i=1:size(tx3,1)
% figure;plot(tx3(i,:))
% [x6(i) y6(i)] = ginput(1);
% close all;
% end
% %% 
% x(x<200)=[];
% x5(x5<200)=[];
% x6(x6<200)=[];
% 
% g1=[];g2=[];g3=[];
% p1=[];p2=[];p3=[];
% p1=x-200;p2=x5-200;p3=x6-200;
% p1(isnan(p1))=[];
% p2(isnan(p2))=[];
% p3(isnan(p3))=[];
% color_id={[1 0 0],([255 165 0]/256),([0 191 255]/256)};
% par=[];par=[p1 p2 p3]';
% g1=[];g1=ones(1,length(p1));
% g2=[];g2=ones(1,length(p2))*2;
% g3=[];g3=ones(1,length(p3))*3;
% gro=[];gro=[g1 g2 g3]';
% fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 250]);set(gcf,'color','w');
% violins = violinplot(par, gro,'ViolinColor',[[0.5 0.5 0.5];[0.5 0.5 0.5];[0.5 0.5 0.5]],'ShowMean', false,'ShowMedian', true,'ViolinAlpha',0,'MedianColor',[0 0 0],'BoxColor',[0 0 0]...
% ,'MarkerSize',10,'MedianMarkerSize',20);
% box off;
% xlim([0 4]);ylabel('First peak time (ms)');set(gca,'FontSize',12);set(gca,'TickDir','out');
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
% text(0.8,-20,'L4','Color','k','FontSize',11);
% text(1.8,-20,'L5','Color','k','FontSize',11);
% text(2.8,-20,'L6','Color','k','FontSize',11);
% 
% %% 
% 
% sdf_tm=[];
% sdf_tm=temp3(:,1:700)-nanmean(temp3(:,1:199),2);
% ani=4;
% color1='m';color2='m';color3='g';color4='r';
% stim_on=20;
% 
% tx1=[];tx2=[];tx3=[];tx4=[];
% tx1=sdf_tm(l5_r{ani},:);
% tx2=sdf_tm(l5_r{ani},:);
% tx3=sdf_tm(l6_r{ani},:);
% tx4=sdf_tm(l4_r{ani},:);
% for i=1:size(tx4,1)
% figure;plot(tx4(i,:))
% [xr(i) yr(i)] = ginput(1);
% close all;
% end
% %% L6
% for i=1:size(tx3,1)
% figure;plot(tx3(i,:))
% [x6r(i) y6r(i)] = ginput(1);
% close all;
% end
% %% 
% %% L5
% for i=1:size(tx1,1)
% figure;plot(tx1(i,:))
% [x5r(i) y5r(i)] = ginput(1);
% close all;
% end
% %% 
% fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 900, 1000]);set(gcf,'color','w');
% t=tiledlayout("vertical",'TileSpacing','Compact');
% nexttile
% plot(tx4(:,100:400)','Color','k');box off
% hold on;line([100 100],[-1 max(get(gca, 'ylim'))],'Color','r');box off;
% title('L4');offsetAxes;
% nexttile
% plot(tx1(:,100:400)','Color','k');box off
% hold on;line([100 100],[-1 max(get(gca, 'ylim'))],'Color','r');box off;
% title('L5');offsetAxes;ylabel('FR (Hz)')
% nexttile
% plot(tx3(:,100:400)','Color','k');box off
% hold on;line([100 100],[-1 max(get(gca, 'ylim'))],'Color','r');box off;
% title('L6');xlabel('Time (ms)');
% offsetAxes;
% %% 
% xr(xr<200)=[];
% x5r(x5r<200)=[];
% x6r(x6r<200)=[];
% g1=[];g2=[];g3=[];
% p1=[];p2=[];p3=[];
% p1=xr-200;p2=x5r-200;p3=x6r-200;
% p1(isnan(p1))=[];
% p2(isnan(p2))=[];
% p3(isnan(p3))=[];
% color_id={[1 0 0],([255 165 0]/256),([0 191 255]/256)};
% par=[];par=[p1 p2 p3]';
% g1=[];g1=ones(1,length(p1));
% g2=[];g2=ones(1,length(p2))*2;
% g3=[];g3=ones(1,length(p3))*3;
% gro=[];gro=[g1 g2 g3]';
% fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 250]);set(gcf,'color','w');
% violins = violinplot(par, gro,'ViolinColor',[[0.5 0.5 0.5];[0.5 0.5 0.5];[0.5 0.5 0.5]],'ShowMean', false,'ShowMedian', true,'ViolinAlpha',0,'MedianColor',[0 0 0],'BoxColor',[0 0 0]...
% ,'MarkerSize',10,'MedianMarkerSize',20);
% box off;
% xlim([0 4]);ylabel('First peak time (ms)');set(gca,'FontSize',12);set(gca,'TickDir','out');
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
% text(0.8,-20,'L4','Color','k','FontSize',11);
% text(1.8,-20,'L5','Color','k','FontSize',11);
% text(2.8,-20,'L6','Color','k','FontSize',11);