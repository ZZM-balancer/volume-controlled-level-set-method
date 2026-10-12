% Ve = 0.0;
% for i=2:Num
%     Ve = Ve + abs(Volume(i)-Volume(1));
% end
% Ve = Ve / Num
figure('Position',[100,100,900,650]);
% 
% plot(0:1:Num-1, Volume(1:Num), "LineWidth", 2.5)
% xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylabel("Volume","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylim([0,1]);

plot(0:1:Num-1, drag(1:Num), "LineWidth", 2.5)
xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylabel("Drag","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylim([10,46]);
% 
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.25;
ax.Box = 'on';
% exportgraphics(gcf,'M2E1_Volume.pdf','ContentType','vector');
exportgraphics(gcf,'M2E1_Drag.pdf','ContentType','vector');