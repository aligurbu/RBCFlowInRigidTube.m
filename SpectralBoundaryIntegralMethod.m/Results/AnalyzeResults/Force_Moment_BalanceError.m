%% Forces and moments balance errors
figure('Color','white')
yyaxis left
plot(T_step,ErrMemForce,'LineWidth',LineWidthind)
xlabel('Time (s)'); ylabel('Force balance error ($\%$)');

yyaxis right
plot(T_step,ErrMemMomentForce,'LineWidth',LineWidthind)
xlabel('Time (s)'); ylabel('Moment balance error ($\%$)')

xlim([T_step(1) T_step(end)]);
legend('Forces error','Moments error','Location','northeast')
legend('boxoff')
set(gca,'FontName','cambria math','FontSize',12)
box on

if verbose_Plot
    set(gcf,'PaperPositionMode','auto')
    print('Force_Moment_BalanceError','-dpng','-r0')
end