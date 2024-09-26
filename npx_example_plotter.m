
function npx_example_plotter(ex_unit1,raster_traces,sdf_all,t_pre,t_post,color_s)

trial_nr=size(raster_traces,1); 

data={};data=raster_traces;

fig2=figure;set(gcf,'color','w');set(fig2, 'Position', [500, 800 ,700, 350]);t=tiledlayout("vertical");
t.TileSpacing = 'compact';t.Padding = 'compact';
nexttile;
ylim([0 10]);
lims=get(gca, 'ylim');

r=rectangle('Position',[0,0,0.5,lims(2)-2],'FaceColor',[color_s 0.8],'EdgeColor',[1 1 1 0]);hold on;
r.FaceAlpha=0.2;
line([0, 0], get(gca, 'ylim')-2, 'color', 'k', 'linestyle', '--','LineWidth',1.5);xlim([-t_pre t_post]);hold on;
for up_i = 1 : trial_nr
s=scatter(data{up_i,ex_unit1}, up_i + ones(length(data{up_i,ex_unit1}), 1));s.Marker='.';s.MarkerEdgeColor='k';
hold on;
end

ylabel('Trials');set(gca,'FontSize',12);
yticks([1:5:11])
box off;h = gca;h.XAxis.Visible = 'off';set(gca,'TickDir','out'); 
  line([max(get(gca, 'xlim'))-0.05 max(get(gca, 'xlim'))],[12 12],'Color','k');
text(max(get(gca, 'xlim'))-0.05,15,'50ms');


data_sdf=[];data_sdf=sdf_all;
tstep=0.001;
time=tstep-t_pre:tstep:t_post;
yticklabels({'1','5','10'});set(gca,'FontSize',12);


nexttile;
plot(time,data_sdf(ex_unit1,[1:round(((t_pre+t_post)*1000))]),'Color',color_s,'LineWidth',1.5);box off;
line([0, 0], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--','LineWidth',1.5);xlim([-t_pre t_post]);
ylabel([]);set(gca,'FontSize',12);
xlabel([]);
%ylim([-15 max(get(gca,'ylim'))])
set(gca,'TickDir','out'); 
xticklabels({[],[],[],[]});
 h = gca;h.YAxis.Visible = 'off';  
  h = gca;h.XAxis.Visible = 'off'; 

line([max(get(gca, 'xlim')) max(get(gca, 'xlim'))],[-5 15],'Color','k');
tt=text(max(get(gca, 'xlim'))+0.03,max(get(gca,'ylim'))-max(get(gca,'ylim'))/2,'20Hz');
 set(tt,'Rotation',270);
end