%% Three axes of ellipsoid with mass and volume of RBC
figure('Color','white')
hold on
plot(T_step, dA, 'k-', T_step, dB, 'b-', T_step, dC, 'g-', 'LineWidth', LineWidthind),
xlabel('Time (s)'), ylabel('Three axes of ellipsoid ($\mu$m)'),
xlim([T_step(1) T_step(end)]);
legend('a','b','c','Location','southeast')
legend('boxoff')
set(gca,'FontName','cambria math','FontSize',12)
box on

if verbose_Plot
    set(gcf,'PaperPositionMode','auto')
    print('AxesOfEllipsoid','-dpng','-r0')
end