function plot_functions(data_cell,type,color,style,lim_x,data_cell2)

    if isequal(type,'PNS')
        plot_norm_star(data_cell,color,style,lim_x);
    elseif isequal(type,'PRS')
        plot_rho_star(data_cell,data_cell2,color,style,lim_x);  
    end


end


function plot_norm_star(data_cell,color,style,lim_x)

    x_int=[1:lim_x];
    q_5=quantile(cell2mat(data_cell'),0.5);
    q_25=quantile(cell2mat(data_cell'),0.25);
    q_75=quantile(cell2mat(data_cell'),0.75);
    loglog(q_5,'Color',color,'Linestyle',style,'LineWidth',2);

    hold on
    fill([x_int,fliplr(x_int)],[q_5,fliplr(q_75)],color,'FaceAlpha','0.2','LineStyle','none')
    fill([x_int,fliplr(x_int)],[q_5,fliplr(q_25)],color,'FaceAlpha','0.2','LineStyle','none')
    xlim([1,inf])
    ylim([1e-10,inf])

end

function plot_rho_star(data_cell,data_cell2,color,style,lim_x)

    x_int=[1:lim_x];
    q_5=quantile(cell2mat(data_cell2')-cell2mat(data_cell'),0.5);
    q_25=quantile(cell2mat(data_cell2')-cell2mat(data_cell'),0.25);
    q_75=quantile(cell2mat(data_cell2')-cell2mat(data_cell'),0.75);
    loglog(q_5,strcat(color,style),'LineWidth',2);

    hold on
    fill([x_int,fliplr(x_int)],[q_5(x_int),fliplr(q_75(x_int))],color,'FaceAlpha','0.2','LineStyle','none')
    fill([x_int,fliplr(x_int)],[q_5(x_int),fliplr(q_25(x_int))],color,'FaceAlpha','0.2','LineStyle','none')
    xlim([1,inf])
    ylim([1e-10,inf])

end
