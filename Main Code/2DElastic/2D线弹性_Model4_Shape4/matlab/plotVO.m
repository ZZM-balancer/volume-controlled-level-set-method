Num = 1001;
figure('Position',[100,100,900,650]);

% plot(0:1:Num-1, Volume(1:Num), "LineWidth", 2.5)
% xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylabel("Volume $\int_D M_{\gamma}(\phi) dx$","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylim([0,2]);


plot(0:1:Num-1, compliance(1:Num), "LineWidth", 2.5)
xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylabel("Compliance $J(\phi)$","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylim([0,15]);
% 
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.25;
ax.Box = 'on';
% exportgraphics(gcf,'M1E4_volume.pdf','ContentType','vector');
exportgraphics(gcf,'M1E4_compliance.pdf','ContentType','vector');
