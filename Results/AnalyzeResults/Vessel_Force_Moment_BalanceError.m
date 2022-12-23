%% Vessel Forces and moments balance errors
figure('Color','white')
yyaxis left
plot(T_step,ForceBalanceError,'LineWidth',LineWidthind)
xlabel('Time (s)'); ylabel('Vessel force balance error ($\%$)');

yyaxis right
plot(T_step,MomentBalanceError,'LineWidth',LineWidthind)
xlabel('Time (s)'); ylabel('Vessel moment balance error ($\%$)')

legend('Forces error','Moments error','Location','southeast')
legend('boxoff')
set(gca,'FontName','cambria math','FontSize',12)
box on

if verbose_Plot
    set(gcf,'PaperPositionMode','auto')
    print('Vessel_Force_Moment_BalanceError','-dpng','-r0')
end