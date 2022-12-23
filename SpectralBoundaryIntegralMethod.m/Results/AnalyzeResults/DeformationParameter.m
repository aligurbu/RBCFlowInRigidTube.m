%% Deformation index (Taylor deformation parameter)
figure('Color','white')
plot(T_step,DeformationIndex,'LineWidth',LineWidthind)
xlabel('Time (s)')
ylabel('Deformation parameter')
set(gca,'FontName','cambria math','FontSize',12)
xlim([T_step(1) T_step(end)]);
box on

if verbose_Plot
    set(gcf,'PaperPositionMode','auto')
    print('DeformationParameter','-dpng','-r0')
end