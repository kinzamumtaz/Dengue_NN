for i=1:17
    initial(i)=1e5;
end

initial(1)=1e6;
initial(15)=1e6;
%for i=2:9
 %   initial(i)=2e4;
%end

[t,y]=ode45(@Dengue3, [0 350], initial);


Y1=log(y(:,3)+y(:,13));
Y2=log(y(:,7)+y(:,9));

Y3=log(y(:,3)+y(:,13)+y(:,7)+y(:,9));

for i=1:17
    min(y(:,i));
end
set(gcf, 'Position', [400 400 600 300]); % e.g., 600x400 pixels

plot(t,Y1, 'r')
hold on;
plot(t,Y2, 'b')

xlabel('time / days','fontsize',10);
ylabel('Log(I_i+ I_{ji})','fontsize',10);
title('R_2= 0.5905 < R_1= 0.6602 < 1', 'fontsize', 10)
legend('Log(I_1 + I_{21})', 'Log(I_2 + I_{12})')
box on;
set(gca, 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k') 
%axes('fontsize', 20)