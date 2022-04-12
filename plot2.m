function plot2(fitness1,fitness2,scaled1,scaled2,iGeneration) 
    f = figure('visible', 'off');
    subplot(2,1,1);
    plot(sort(fitness1,'descend'),'ro');hold on;
    plot(sort(fitness2,'descend'),'go');
    subplot(2,1,2);
    plot(sort(scaled1,'descend'),'ro');hold on;
    plot(sort(scaled2,'descend'),'go');
    titlename = ['iteration time=',num2str(iGeneration),'.png'];
    saveas(f,titlename);
end