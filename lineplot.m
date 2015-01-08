function lineplot(par,clr)

m = - par(1, :)./par(2, :);
c = - par(3, :)./par(2, :);
% th = 0.1; 

% figure(h); 
% subplot(2,3,2); 
hold on; 
for i=1:numel(m)
    x = linspace(-1,1,100);
    y = m(i)*x + repmat(c(i),1,100);
    idx = abs(y) < 1; 
    plot(x(idx), y(idx),'-','color',clr(i+1, :)); 
    y = linspace(-1,1,100);
    x = y/m(i) - repmat(c(i)/m(i),1,100);
    idx = abs(x) < 1; 
    plot(x(idx), y(idx),'-','color',clr(i+1, :)); 
    
    
end
hold off; 

% y = m*x + repmat(c-th,1,100);
% plot(x,y,'--','color',clr);
% 
% y = m*x + repmat(c+th,1,100);
% plot(x,y,'--','color',clr);
end