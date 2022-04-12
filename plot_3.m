function plot_3(fitness1,fitness2,parents1,parents2, pselected1,pselected2,iGeneration) 
%     f = figure('visible', 'off');
    subplot(3,1,1);
    plot(sort(fitness1,'descend'),'ro');hold on;
    plot(sort(fitness2,'descend'),'go');
    subplot(3,1,2);
    plot(sort(parents1,'descend'),'ro');hold on;
    plot(sort(parents2,'descend'),'go');
    subplot(3,1,3);
    plot(sort(pselected1,'descend'),'ro');hold on;
    plot(sort(pselected2,'descend'),'go');
%     titlename = ['iteration time=',num2str(iGeneration),'.png'];
%     saveas(f,titlename);
end