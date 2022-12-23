%% Area Error and Volume Error
figure('Color','white')
yyaxis left
plot(T_step,Area,'LineWidth',LineWidthind), 
xlabel('Time (s)'), ylabel('Area ($\mu m^2$)'),

yyaxis right
plot(T_step,Volume,'LineWidth',LineWidthind), 
xlabel('Time (s)'), ylabel('Volume ($\mu m^3$)'),
% title(['Area Error = ', num2str(ErrorArea),' $\%$',...
%        '       Volume Error = ', num2str(ErrorVolume),' $\%$'])
xlim([T_step(1) T_step(end)]);
legend('Area','Volume','Location','east')
legend('boxoff')
set(gca,'FontName','cambria math','FontSize',12)
box on

if verbose_Plot
    set(gcf,'PaperPositionMode','auto')
    print('Area_Volume_Change','-dpng','-r0')
end