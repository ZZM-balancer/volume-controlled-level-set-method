Num = 701;
figure('Position',[100,100,900,650]);

plot(0:1:Num-1, Volume(1:Num), "LineWidth", 2.5)
xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylabel("Volume $\int_D M_{\gamma}(\phi) dx$","Interpreter","latex","FontSize",18,'FontWeight','bold');
ylim([0,2.5]);

% plot(0:1:Num-1, lambda1(1:Num), "LineWidth", 2.5)
% xlabel("Iteration number $n$","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylabel("$\lambda_1$","Interpreter","latex","FontSize",18,'FontWeight','bold');
% ylim([0,2]);

ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.25;
ax.Box = 'on';
exportgraphics(gcf,'M5E2_Volume.pdf','ContentType','vector');
% exportgraphics(gcf,'M5E2_lambda.pdf','ContentType','vector');