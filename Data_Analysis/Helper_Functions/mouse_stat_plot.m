function mouse_stat_plot(temp_mat)

conds = size(temp_mat,1);
mice = size(temp_mat,2);
temp_v_1 = repmat([1:conds]',1,mice);
temp_v_2 = repmat(linspace(-0.15,0.15,mice),conds,1);
plotvect = temp_v_1+temp_v_2;
plotvect = plotvect(:);

median_forbar = nanmean(temp_mat');
bar([1:conds],median_forbar,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
hold on
scatter(plotvect,temp_mat(:),20,round(plotvect),'filled')
caxis([0 conds+.5])
xlim([0.5 conds+.5])
ylim([min(0,min(temp_mat(:))*1.1) max(0,max(temp_mat(:))*1.1+0.001)])
colormap('jet')