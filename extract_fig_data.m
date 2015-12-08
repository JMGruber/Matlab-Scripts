h=gcf
axesObjs = get(h, 'Children'); %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata1 = get(dataObjs{2}(1), 'XData');
ydata1 = get(dataObjs{2}(1), 'YData');
xdata2 = get(dataObjs{2}(2), 'XData');
ydata2 = get(dataObjs{2}(2), 'YData');
ydata3 = get(dataObjs{2}(3), 'YData');
xdata3 = get(dataObjs{2}(3), 'XData');
figure(45)
plot(xdata1,ydata1)
hold on
plot(xdata2,ydata2)
plot(xdata3,ydata3)