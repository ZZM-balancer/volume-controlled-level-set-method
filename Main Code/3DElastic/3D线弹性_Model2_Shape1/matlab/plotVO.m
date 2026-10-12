% plot(0:1:Num-1, Volume(1:Num), "LineWidth", 1.5)
% xlabel("Iteration number $n$","Interpreter","latex");
% ylabel("Volume $\int_D H_{\gamma}(\phi) dx$","Interpreter","latex");
% ylim([0,1]);


plot(0:1:Num-1, compliance(1:Num), "LineWidth", 1.5)
xlabel("Iteration number $n$","Interpreter","latex");
ylabel("Compliance $J(\Omega)$","Interpreter","latex");
ylim([0,100]);