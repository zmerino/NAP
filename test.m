
% 
qnum = 20

% disp(sprintf('$Q_%s(x)$',[num2str(100/qnum),'%']),'Interpreter','latex')

x = linspace(0,1,10);


figure()
plot(x,x,'-r', 'DisplayName',sprintf('$Q_{%s}(x)$',[num2str(100/qnum),'\%']))
legend('Interpreter','latex')


disp(sprintf("$Q_%s(x)$",num2str(100/qnum)))