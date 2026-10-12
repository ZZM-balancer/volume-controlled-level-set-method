figure('Position',[100,100,900,650]);
% Ve = 0.0;
% for i=2:Num
%     Ve = Ve + abs(Volume(i)-Volume(1));
% end
% Ve = Ve / Num

% plot(0:1:Num-1, Volume(1:Num), "LineWidth", 2.5)
% xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylabel("Volume","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylim([0.6,1.5]);


plot(0:1:Num-1, compliance(1:Num), "LineWidth", 2.5)
xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylabel("Compliance","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylim([0,4]);
% 
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.25;
ax.Box = 'on';
% exportgraphics(gcf,'M1E1_volume.pdf','ContentType','vector');
exportgraphics(gcf,'M1E1_compliance.pdf','ContentType','vector');