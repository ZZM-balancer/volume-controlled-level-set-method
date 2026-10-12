% Volume1 = Volume; compliance1 = compliance;
% Volume2 = Volume; compliance2 = compliance;
% Ve = 0.0;
% for i=2:Num
%     Ve = Ve + abs(Volume(i)-Volume(1));
% end
% Ve = Ve / Num
%% ===================== Volume =====================
% figure('Position',[100,100,900,650]);
% 
% plot(0:Num-1, Volume1(1:Num), 'b-', 'LineWidth', 2.5);
% hold on;
% plot(0:Num-1, Volume2(1:Num), '--', 'Color', [0.5, 0, 0.5], 'LineWidth', 2.5);
% hold off;
% 
% xlabel('Iteration number $n$', 'Interpreter', 'latex', ...
%     'FontSize', 18, 'FontWeight', 'bold');
% ylabel('Volume', 'Interpreter', 'latex', ...
%     'FontSize', 18, 'FontWeight', 'bold');
% legend('first initial shape', 'second initial shape', 'Location', 'best', ...
%     'FontSize', 16);
% ylim([0,1]);
% 
% ax = gca;
% ax.FontSize = 18;
% ax.LineWidth = 1.25;
% ax.Box = 'on';
% exportgraphics(gcf,'M4E2_volume.pdf','ContentType','vector');


%% ===================== Compliance =====================
figure('Position',[100,100,900,650]);

plot(0:Num-1, compliance1(1:Num), 'LineWidth', 2.5);
hold on;
plot(0:Num-1, compliance2(1:Num), '--', 'Color', [0.5, 0, 0.5], 'LineWidth', 2.5);
hold off;

xlabel('Iteration number $n$', 'Interpreter', 'latex', ...
    'FontSize', 18, 'FontWeight', 'bold');
ylabel('Compliance', 'Interpreter', 'latex', ...
    'FontSize', 18, 'FontWeight', 'bold');
legend('first initial shape', 'second initial shape', 'Location', 'best', ...
    'FontSize', 16);
ylim([0,100]);

ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.25;
ax.Box = 'on';
exportgraphics(gcf,'M4E2_compliance.pdf','ContentType','vector');