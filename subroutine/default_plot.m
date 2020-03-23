function default_plot(bool_log)

if bool_log
    set(gca,'XScale','log','YScale','log','MinorGridLineStyle','none','FontSize',16);
    grid on
else
    set(gca,'MinorGridLineStyle','none','FontSize',16);
    grid on
end
