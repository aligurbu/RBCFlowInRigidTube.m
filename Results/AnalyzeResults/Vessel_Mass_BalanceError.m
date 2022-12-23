%% Vessel mass balance and mass balance error
%% MassBalance(nframe) = (inflow-outflow);
%% MassBalanceError(nframe) = ((inflow-outflow)/(inflow+outflow))*100;
figure('Color','white')
yyaxis left
plot(T_step,TotalInflow, ...
     T_step,TotalDischarge,'k','LineWidth',LineWidthind)
xlabel('Time (s)'); ylabel('Total inflow and outflow ($\mu$m$^3$/s)');

yyaxis right
plot(T_step,MassBalanceError,'LineWidth',LineWidthind)
xlabel('Time (s)'); ylabel('Vessel mass balance error ($\%$)')

xlim([T_step(1) T_step(end)]);
legend('Inflow','Outflow','Mass error','Location','southeast')
legend('boxoff')
set(gca,'FontName','cambria math','FontSize',12)
box on

if verbose_Plot
    set(gcf,'PaperPositionMode','auto')
    print('Vessel_Mass_Balance_Error','-dpng','-r0')
end