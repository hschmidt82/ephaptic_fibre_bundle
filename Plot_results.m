% This script loads data from Results.mat and plots them as in Figure 8 in
% Schmidt, Hahn, Deco, Knoesche; PLOS Computational Biology (2020).

close all; clear all;

load('Results.mat')

dmat = cat(5,dmat0,dmat);
dmat2 = cat(5,dmat02,dmat2);
lmat = cat(5,lmat0,lmat);

x = [0.1:0.1:1];
delay2 = squeeze(mean(dmat,2));
delay2std = squeeze(std(dmat,[],2));
t3max2 = squeeze(mean(lmat,2));
t3max2std = squeeze(std(lmat,[],2));
delay4 = squeeze(mean(dmat2,2));
delay4std = squeeze(std(dmat2,[],2));

c1(1,:) = [0 84 95]./255;
c1(2,:) = [163 2 52]./255;
c1(3,:) = [227 124 29]./255;
c1(4,:) = [103 119 26]./255;
alpha_val = 0.14;

for o = 1:3

% panel 1:
figure;
for m = 1:4
    
y1 = squeeze(delay2(o,:,m,2)+delay2std(o,:,m,2));
y2 = squeeze(delay2(o,:,m,2)-delay2std(o,:,m,2));
y1 = y1./100; y2 = y2./100;

s = fill([x x(end:-1:1)],[y1 y2(end:-1:1)],c1(m,:),'EdgeColor','None');
alpha(s,alpha_val);
hold on;
plot(x,0.5.*(y1+y2),'Color',c1(m,:),'LineWidth',2)

end

for m = 1:4
    
y1 = squeeze(delay2(o,:,m,1)+delay2std(o,:,m,1));
y2 = squeeze(delay2(o,:,m,1)-delay2std(o,:,m,1));
y1 = y1./100; y2 = y2./100;

s = fill([x x(end:-1:1)],[y1 y2(end:-1:1)],c1(m,:),'EdgeColor','None');
alpha(s,alpha_val);
hold on;
plot(x,0.5.*(y1+y2),'--','Color',c1(m,:),'LineWidth',2)

end

xlim([0.1 1])
ylim([20 38])

% panel 2:
figure;
for m = 1:4
    
y1 = squeeze(t3max2(o,:,m,2)+t3max2std(o,:,m,2));
y2 = squeeze(t3max2(o,:,m,2)-t3max2std(o,:,m,2));
y1 = y1./100; y2 = y2./100;

s = fill([x x(end:-1:1)],[y1 y2(end:-1:1)],c1(m,:),'EdgeColor','None');
alpha(s,alpha_val);
hold on;
plot(x,0.5.*(y1+y2),'Color',c1(m,:),'LineWidth',2)

end

for m = 1:4
    
y1 = squeeze(t3max2(o,:,m,1)+t3max2std(o,:,m,1));
y2 = squeeze(t3max2(o,:,m,1)-t3max2std(o,:,m,1));
y1 = y1./100; y2 = y2./100;

s = fill([x x(end:-1:1)],[y1 y2(end:-1:1)],c1(m,:),'EdgeColor','None');
alpha(s,alpha_val);
hold on;
plot(x,0.5.*(y1+y2),'--','Color',c1(m,:),'LineWidth',2)

end

xlim([0.1 1])
ylim([36 50])

% panel 3:
figure;
for m = 1:4
    
y1 = squeeze(delay4(o,:,m,2)+delay4std(o,:,m,2));
y2 = squeeze(delay4(o,:,m,2)-delay4std(o,:,m,2));
y1 = y1./100; y2 = y2./100;

s = fill([x x(end:-1:1)],[y1 y2(end:-1:1)],c1(m,:),'EdgeColor','None');
alpha(s,alpha_val);
hold on;
plot(x,0.5.*(y1+y2),'Color',c1(m,:),'LineWidth',2)

end

for m = 1:4
    
y1 = squeeze(delay4(o,:,m,1)+delay4std(o,:,m,1));
y2 = squeeze(delay4(o,:,m,1)-delay4std(o,:,m,1));
y1 = y1./100; y2 = y2./100;

s = fill([x x(end:-1:1)],[y1 y2(end:-1:1)],c1(m,:),'EdgeColor','None');
alpha(s,alpha_val);
hold on;
plot(x,0.5.*(y1+y2),'--','Color',c1(m,:),'LineWidth',2)

end

xlim([0.1 1])
ylim([13 19])

end

% END OF SCRIPT