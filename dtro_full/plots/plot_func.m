clear all
close all
%%
% Plot convergence in objective for different trees

data_topo=load('../topo.mat');
%data_topo_er=load('topo_er.mat');
%data_topo.norm_star_cell_ER_03=data_topo_er.norm_star_cell_ER_03;
%data_topo.norm_star_cell_ER_06=data_topo_er.norm_star_cell_ER_06;
%data_topo.norm_star_cell_ER_09=data_topo_er.norm_star_cell_ER_09;


plot_functions(data_topo.norm_star_cell_TP,'PNS','r','-',length(data_topo.norm_star_cell_TP{1}));
plot_functions(data_topo.norm_star_cell_TR2,'PNS','c','-',length(data_topo.norm_star_cell_TR2{1}));
plot_functions(data_topo.norm_star_cell_RT,'PNS',[0.4660, 0.6740, 0.1880],'-',length(data_topo.norm_star_cell_RT{1}));
plot_functions(data_topo.norm_star_cell_TR5,'PNS','m','-',length(data_topo.norm_star_cell_TR5{1}));
plot_functions(data_topo.norm_star_cell_TS,'PNS','b','-',length(data_topo.norm_star_cell_TS{1}));
plot_functions(data_topo.norm_star_cell_ER_03,'PNS','g','-',length(data_topo.norm_star_cell_ER_03{1}));
plot_functions(data_topo.norm_star_cell_ER_06,'PNS',[0.8500, 0.3250, 0.0980],'-',length(data_topo.norm_star_cell_ER_06{1}));
plot_functions(data_topo.norm_star_cell_ER_09,'PNS',[0.8500, 0.6, 0.0980],'-',length(data_topo.norm_star_cell_ER_09{1}));
plot_functions(data_topo.norm_star_cell_FC,'PNS','k','-',length(data_topo.norm_star_cell_FC{1}));
grid on
ylim([1e-6,inf])


p1=loglog(rand,'r-','LineWidth',2,'Visible','off');
p2=loglog(rand,'c-','LineWidth',2,'Visible','off');
p3=loglog(rand,'Color',[0.4660, 0.6740, 0.1880],'Linestyle','-','LineWidth',2,'Visible','off');
p4=loglog(rand,'m-','LineWidth',2,'Visible','off');
p5=loglog(rand,'b-','LineWidth',2,'Visible','off');
p6=loglog(rand,'Color','g','Linestyle','-','LineWidth',2,'Visible','off');
p7=loglog(rand,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','-','LineWidth',2,'Visible','off');
p8=loglog(rand,'Color',[0.8500, 0.6, 0.0980],'Linestyle','-','LineWidth',2,'Visible','off');
p9=loglog(rand,'k-','LineWidth',2,'Visible','off');

nbnodes=31;
addpath('/Users/musluoglucem/Documents/PhD/DTRO/d_trace_ratio/dtro_full')
ac_tp=alg_connect('TP',31);
ac_tr2=alg_connect('TR',31,2);
ac_rt=median(cell2mat(data_topo.alg_c_RT));
ac_tr5=alg_connect('TR',31,5);
ac_ts=alg_connect('TS',31);
ac_er_03=median(cell2mat(data_topo.alg_c_ER_03));
ac_er_06=median(cell2mat(data_topo.alg_c_ER_06));
ac_er_09=median(cell2mat(data_topo.alg_c_ER_09));
ac_fc=alg_connect('FC',31);

legend_tp=[char("Path,\hskip 2em AC: ") num2str(ac_tp,'%.3f')];
legend_tr2=[char("Reg. 2,\hskip 1.25em AC: ") num2str(ac_tr2,'%.3f')];
legend_rt=[char("Rand.,\hskip 1.5em AC: ") num2str(ac_rt,'%.3f')];
legend_tr5=[char("Reg. 5,\hskip 1.25em AC: ") num2str(ac_tr5,'%.3f')];
legend_ts=[char("Star,\hskip 2.25em AC: ") num2str(ac_ts,'%.3f')];
legend_er_03=[char("ER(0.3),\hskip 0.7em AC: ") num2str(ac_er_03,'%.3f')];
legend_er_06=[char("ER(0.6),\hskip 0.7em AC: ") num2str(ac_er_06,'%.3f')];
legend_er_09=[char("ER(0.9),\hskip 0.7em AC: ") num2str(ac_er_09,'%.3f')];
legend_fc=[char("FC,\hskip 2.75em AC: ") num2str(ac_fc,'%.3f')];

[hl(1).leg, hl(1).obj, hl(1).hout, hl(1).mout] = ...
legendflex([p1;p2;p3;p4;p5;p6;p7;p8;p9], {legend_tp,legend_tr2,legend_rt,legend_tr5,...
legend_ts,legend_er_03,legend_er_06,legend_er_09,legend_fc}, ...
'anchor', {'sw','sw'}, ...
'buffer', [7 7], ...
'fontsize',12, ...
'title', 'Network Type',...
'interpreter','latex');
%xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
ylabel('$\epsilon(X^i)$','Interpreter','latex','FontSize', 20)
xlabel('Iterations $i$','Interpreter','latex','FontSize', 20)
set(gcf,'Position',[440   378   560   380])
%%

qeq7=load('../Q7.mat');
qeq5=load('../Q5.mat');
qeq3=load('../Q3.mat');
qeq1=load('../Q1.mat');

% Plot convergence in argument

plot_functions(qeq7.norm_star_cell_FC,'PNS','c','-',length(qeq7.norm_star_cell_FC{1}));
plot_functions(qeq7.norm_star_cell_RT,'PNS','c','--',length(qeq7.norm_star_cell_RT{1}));
plot_functions(qeq5.norm_star_cell_FC,'PNS','g','-',length(qeq5.norm_star_cell_FC{1}));
plot_functions(qeq5.norm_star_cell_RT,'PNS','g','--',length(qeq5.norm_star_cell_RT{1}));
plot_functions(qeq3.norm_star_cell_FC,'PNS','b','-',length(qeq3.norm_star_cell_FC{1}));
plot_functions(qeq3.norm_star_cell_RT,'PNS','b','--',length(qeq3.norm_star_cell_RT{1}));
plot_functions(qeq1.norm_star_cell_FC,'PNS','r','-',length(qeq1.norm_star_cell_FC{1}));
plot_functions(qeq1.norm_star_cell_RT,'PNS','r','--',length(qeq1.norm_star_cell_RT{1}));
grid on

p1=loglog(rand,'k-','LineWidth',2,'Visible','off');
p2=loglog(rand,'k--','LineWidth',2,'Visible','off');
p3=fill([1 2],[1 2],'c','LineStyle','none','Visible','off');
p4=fill([1 2],[1 2],'g','LineStyle','none','Visible','off');
p5=fill([1 2],[1 2],'b','LineStyle','none','Visible','off');
p6=fill([1 2],[1 2],'r','LineStyle','none','Visible','off');

[hl(1).leg, hl(1).obj, hl(1).hout, hl(1).mout] = ...
legendflex([p1;p2], {'Fully Connected','Random Tree'}, ...
'anchor', {'sw','sw'}, ...
'buffer', [7 7], ...
'fontsize',12, ...
'title', 'Network Type',...
'interpreter','latex');
[hl(2).leg, hl(2).obj, hl(2).hout, hl(2).mout] = ...
legendflex([p6;p5;p4;p3], {'$Q=1$','$Q=3$','$Q=5$','$Q=7$'}, ...
'ref', hl(1).leg, ...
'anchor', {'se','sw'}, ...
'ncol', 2,...
'buffer', [0 0], ...
'fontsize', 12', ...
'title', 'Projection Dimension',...
'interpreter','latex');
%xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
ylabel('$\epsilon(X^i)$','Interpreter','latex','FontSize', 20)
xlabel('Iterations $i$','Interpreter','latex','FontSize', 20)
set(gcf,'Position',[440   378   560   280])

%%

qeq10=load('../K10.mat');
qeq50=load('../K50.mat');

plot_functions(qeq10.norm_star_cell_FC,'PNS','r','-',length(qeq10.norm_star_cell_FC{1}));
plot_functions(qeq10.norm_star_cell_RT,'PNS','r','--',length(qeq10.norm_star_cell_RT{1}));
plot_functions(qeq3.norm_star_cell_FC,'PNS','b','-',length(qeq3.norm_star_cell_FC{1}));
plot_functions(qeq3.norm_star_cell_RT,'PNS','b','--',length(qeq3.norm_star_cell_RT{1}));
plot_functions(qeq50.norm_star_cell_FC,'PNS','g','-',length(qeq50.norm_star_cell_FC{1}));
plot_functions(qeq50.norm_star_cell_RT,'PNS','g','--',length(qeq50.norm_star_cell_RT{1}));
grid on

p1=loglog(rand,'k-','LineWidth',2,'Visible','off');
p2=loglog(rand,'k--','LineWidth',2,'Visible','off');
p3=fill([1 2],[1 2],'r','LineStyle','none','Visible','off');
p4=fill([1 2],[1 2],'b','LineStyle','none','Visible','off');
p5=fill([1 2],[1 2],'g','LineStyle','none','Visible','off');

[hl(1).leg, hl(1).obj, hl(1).hout, hl(1).mout] = ...
legendflex([p1;p2], {'Fully Connected','Random Tree'}, ...
'anchor', {'sw','sw'}, ...
'buffer', [7 7], ...
'fontsize',12, ...
'title', 'Network Type',...
'interpreter','latex');
[hl(2).leg, hl(2).obj, hl(2).hout, hl(2).mout] = ...
legendflex([p3;p4;p5], {'K=10','K=30','K=50'}, ...
'ref', hl(1).leg, ...
'anchor', {'se','sw'}, ...
'ncol', 2,...
'buffer', [0 0], ...
'fontsize', 12', ...
'title', 'Number of Nodes',...
'interpreter','latex');
%xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
ylabel('$\epsilon(X^i)$','Interpreter','latex','FontSize', 20)
xlabel('Iterations $i$','Interpreter','latex','FontSize', 20)
set(gcf,'Position',[440   378   560   280])

%%

% Plot convergence rate

conv_rate=load('../conv_rate.mat');

x_int=[2:3000];
data_mat=cell2mat(conv_rate.norm_star_cell_RT');
conv_rates=sqrt(abs(data_mat(:,2:end))./abs(data_mat(:,1:end-1)));
q_5=quantile(conv_rates,0.5);
q_25=quantile(conv_rates,0.25);
q_75=quantile(conv_rates,0.75);
figure
p3=loglog([2:3000],q_5,'b-','LineWidth',2);
xlim([1,inf])
%ylim([0.98,1.02])
grid on
hold on
%p4=semilogx([2:3000],q_25,'b--','LineWidth',0.3);
%semilogx([2:3000],q_75,'b--','LineWidth',0.3)
fill([x_int,fliplr(x_int)],[q_5,fliplr(q_25)],'b','FaceAlpha','0.2','LineStyle','none')
fill([x_int,fliplr(x_int)],[q_5,fliplr(q_75)],'b','FaceAlpha','0.2','LineStyle','none')

data_mat=cell2mat(conv_rate.norm_star_cell_FC');
conv_rates=sqrt(abs(data_mat(:,2:end))./abs(data_mat(:,1:end-1)));
q_5=quantile(conv_rates,0.5);
q_25=quantile(conv_rates,0.25);
q_75=quantile(conv_rates,0.75);
p1=loglog([2:3000],q_5,'r-','LineWidth',2);
xlim([1,inf])
%ylim([0.98,1.02])
%p2=semilogx([2:3000],q_25,'r--','LineWidth',0.3);
%semilogx([2:3000],q_75,'r--','LineWidth',0.3)
fill([x_int,fliplr(x_int)],[q_5,fliplr(q_25)],'r','FaceAlpha','0.2','LineStyle','none')
fill([x_int,fliplr(x_int)],[q_5,fliplr(q_75)],'r','FaceAlpha','0.2','LineStyle','none')

%{
[hl(1).leg, hl(1).obj, hl(1).hout, hl(1).mout] = ...
legendflex([p1;p2], {'Median','25th and 75th percentile'}, ...
'anchor', {'s','se'}, ...
'buffer', [80 7], ...
'fontsize',8, ...
'title', 'FC',...
'interpreter','latex');
[hl(2).leg, hl(2).obj, hl(2).hout, hl(2).mout] = ...
legendflex([p3;p4], {'Median', '25th and 75th percentile'}, ...
'ref', hl(1).leg, ...
'anchor', {'se','sw'}, ...
'buffer', [0 0], ...
'fontsize', 8', ...
'title', 'TG',...
'interpreter','latex');
%xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
%}
[hl(1).leg, hl(1).obj, hl(1).hout, hl(1).mout] = ...
legendflex([p1;p3], {'Fully Connected','Random Tree'}, ...
'anchor', {'s','sw'}, ...
'buffer', [80 7], ...
'fontsize',12, ...
'title', 'Network type',...
'interpreter','latex');

ylabel('$r(X^i)$','Interpreter','latex','FontSize',20)
xlabel('Iterations $i$','Interpreter','latex','FontSize', 20)
set(gcf,'Position',[440   378   560   280])



