%%% Relative level
clear all;
close all;

% Sub 1
% EIL
overall_real_eil = [-1.4 6.3 -37.9 4.6];
trees_real_eil = [2.7 15 -52.4 5.4];
hobble_real_eil = [1.9 4.7 -39.8 -1.6];
crocus_real_eil = [-4.1 10.2 -41.9 4.3];
spring_real_eil = [3.7 26 -56.5 15.8];

overall_ls_eil = [0.9 8.1 -42.5 0.4];
trees_ls_eil = [0.6 15.3 -57.0 1.3];
hobble_ls_eil = [-0.3 7.6 -44.0 1.0];
crocus_ls_eil = [-0.4 10.6 -43.6 -2.0];
spring_ls_eil = [0.5 29.4 -65.3 5.5];

figure;
subplot(2,2,1);
y = [trees_real_eil(1) - trees_ls_eil(1) trees_real_eil(2) - trees_ls_eil(2) 0 trees_real_eil(4) - trees_ls_eil(4)];
z = y;
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 1, EIL: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Trees: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,2);
y = [hobble_real_eil(1) - hobble_ls_eil(1) hobble_real_eil(2) - hobble_ls_eil(2) 0 hobble_real_eil(4) - hobble_ls_eil(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b')
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Hobble: B');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,3);
y = [crocus_real_eil(1) - crocus_ls_eil(1) crocus_real_eil(2) - crocus_ls_eil(2) 0 crocus_real_eil(4) - crocus_ls_eil(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b')
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Crocus: 2nd C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,4);
y = [spring_real_eil(1) - spring_ls_eil(1) spring_real_eil(2) - spring_ls_eil(2) 0 spring_real_eil(4) - spring_ls_eil(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Spring: P');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});


% ATR
overall_real_atr = [-1.0 7.5 -35.6 2.1];
trod_real_atr = [1.7 24 -56.3 5.8];
crystal_real_atr = [11.2 33.5 -60.8 9.3];
god_real_atr = [8.4 27 -59 -1.7];
beautiful_real_atr = [-0.3 7 -35.8 -1.7];

overall_ls_atr = [1.1 8.9 -40.8 1.0];
trod_ls_atr = [0.6 23 -59.4 3];
crystal_ls_atr = [3.2 30.1 -64 5.8];
god_ls_atr = [2.8 28.9 -64.5 6];
beautiful_ls_atr = [-0.8 7.6 -41.9 1.1];
figure;
subplot(2,2,1);
y = [trod_real_atr(1) - trod_ls_atr(1) trod_real_atr(2) - trod_ls_atr(2) 0 trod_real_atr(4) - trod_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 1, ATR: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Trod: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,2);
y = [crystal_real_atr(1) - crystal_ls_atr(1) crystal_real_atr(2) - crystal_ls_atr(2) 0 crystal_real_atr(4) - crystal_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Crystal: C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,3);
y = [god_real_atr(1) - god_ls_atr(1) god_real_atr(2) - god_ls_atr(2) 0 god_real_atr(4) - god_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('God: G');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,4);
y = [beautiful_real_atr(1) - beautiful_ls_atr(1) beautiful_real_atr(2) - beautiful_ls_atr(2) 0 beautiful_real_atr(4) - beautiful_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Beautiful: B');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

% WNDL
overall_real_wndl = [0.7 8.3 -37.1 5];
den_real_wndl = [11.7 23.5 -53.7 6.5];
lieben_real_wndl = [3.3 28.1 -53.1 9.6];
hoffet_real_wndl = [2.8 12.8 -48.2 1.8];
keit_real_wndl = [4.6 27.1 -58.8 5.4];

overall_ls_wndl = [0.7 8 -41.9 0.2];
den_ls_wndl = [4.7 17.1 -54.5 4.8];
lieben_ls_wndl = [0.3 10.2 -48.1 2.4];
hoffet_ls_wndl = [-0.3 10.6 -52.9 1.3];
keit_ls_wndl = [4.4 25.4 -60.6 2.9];

figure;
subplot(2,2,1);
y = [den_real_wndl(1) - den_ls_wndl(1) den_real_wndl(2) - den_ls_wndl(2) 0 den_real_wndl(4) - den_ls_wndl(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 1, WNDL: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Den: D');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,2);
y = [lieben_real_wndl(1) - lieben_ls_wndl(1) lieben_real_wndl(2) - lieben_ls_wndl(2) 0 lieben_real_wndl(4) - lieben_ls_wndl(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Lieben: B');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,3);
y = [hoffet_real_wndl(1) - hoffet_ls_wndl(1) hoffet_real_wndl(2) - hoffet_ls_wndl(2) 0 hoffet_real_wndl(4) - hoffet_ls_wndl(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Hoffet: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,4);
y = [keit_real_wndl(1) - keit_ls_wndl(1) keit_real_wndl(2) - keit_ls_wndl(2) 0 keit_real_wndl(4) - keit_ls_wndl(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
     text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Traurigkeit: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

%Overall
figure;
subplot(3,1,1);
y = [overall_real_eil(1) - overall_ls_eil(1) overall_real_eil(2) - overall_ls_eil(2) 0 overall_real_eil(4) - overall_ls_eil(4)];
x = y;
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('A) Subject 1: Overall Source Level Variation','fontweight','bold','fontsize',12);
xlabel('EIL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

subplot(3,1,2);
y = [overall_real_atr(1) - overall_ls_atr(1) overall_real_atr(2) - overall_ls_atr(2) 0 overall_real_atr(4) - overall_ls_atr(4)];
x = [x; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('ATR');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

subplot(3,1,3);
y = [overall_real_wndl(1) - overall_ls_wndl(1) overall_real_wndl(2) - overall_ls_wndl(2) 0 overall_real_wndl(4) - overall_ls_wndl(4)];
x = [x; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('WNDL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

% Overall room only
figure;
subplot(3,1,1);
y = [overall_ls_eil(1) overall_ls_eil(4) overall_ls_eil(2)];
plot(y,'o','MarkerFaceColor','b');
title('Sub. 1 - Level Relative to AC: Room Only','fontweight','bold','fontsize',12);
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
xlabel('EIL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});
subplot(3,1,2);
y = [overall_ls_atr(1) overall_ls_atr(4) overall_ls_atr(2)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
xlabel('ATR');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});
subplot(3,1,3);
y = [overall_ls_wndl(1) overall_ls_wndl(4) overall_ls_wndl(2)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
xlabel('WNDL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});

% Averages
figure;
errorbar(mean(z), std(z), 'o','MarkerFaceColor','b');
hold on;
plot(median(z),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('A) Subject 1: Consonant Level Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Level (dB)');
hline = refline([0 0]);
set(hline,'Color','k');

figure;
errorbar(mean(x), std(x), 'o','MarkerFaceColor','b');
hold on;
plot(median(x),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('A) Subject 1: Overall Level Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Level (dB)');
hline = refline([0 0]);
set(hline,'Color','k');


% Sub2
% EIL
overall_real_eil = [-32.0 1.6 2.8 9.4];
trees_real_eil = [-43.3 3.5 2.0 13.5];
hobble_real_eil = [-46.7 7.1 13.6 17.9];
crocus_real_eil = [-38.9 -2.7 0.2 15];
spring_real_eil = [-53.9 6.9 7.4 34.5];

overall_ls_eil = [-36.8 0.5 -0.1 7.5];
trees_ls_eil = [-46.8 0.5 1.5 14.3];
hobble_ls_eil = [-44.4 2.3 3.2 11.8];
crocus_ls_eil = [-43.7 13 2.4 16];
spring_ls_eil = [-55.3 1.1 3.2 19.2];

figure;
subplot(2,2,1);
y = [0 trees_real_eil(2:end) - trees_ls_eil(2:end)];
z = [y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 2, EIL: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Trees: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,2);
y = [0 hobble_real_eil(2:end) - hobble_ls_eil(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Hobble: B');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,3);
y = [0 crocus_real_eil(2:end) - crocus_ls_eil(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Crocus: 2nd C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,4);
y = [0 spring_real_eil(2:end) - spring_ls_eil(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Spring: P');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});


% ATR
overall_real_atr = [-32.4 1.5 2.8 11.4];
trod_real_atr = [-52.9 2.7 5.9 24.5];
crystal_real_atr = [-49.3 7.8 -0.4 23.5];
god_real_atr = [-53.0 1.4 2.3 23.4];
beautiful_real_atr = [-47.1 8.5 8.0 10.6];

overall_ls_atr = [-36.3 1.3 0.8 7.6];
trod_ls_atr = [-56.5 1 2.6 23.8];
crystal_ls_atr = [-52.8 -0.1 0.9 19.5];
god_ls_atr = [-55.8 0.2 3.5 20.2];
beautiful_ls_atr = [-51.2 0.4 4.4 21.1];

figure;
subplot(2,2,1);
y = [0 trod_real_atr(2:end) - trod_ls_atr(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 2, ATR: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Trod: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,2);
y = [0 crystal_real_atr(2:end) - crystal_ls_atr(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Crystal: C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,3);
y = [0 god_real_atr(2:end) - god_ls_atr(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('God: G');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,4);
y = [0 beautiful_real_atr(2:end) - beautiful_ls_atr(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Beautiful: B');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});


% JTOU
overall_real_jtou = [-31.0 1.6 2.6 8.2];
crystal_real_jtou = [-45.9 5.1 3.6 18.3];
spend_real_jtou = [-45.9 -16.2 6.7 19.2];
two_real_jtou = [-41.0 2.3 -2.2 16.0];
just_real_jtou = [-41.3 0.1 4.2 21.0];

overall_ls_jtou = [-34.8 0.2 0.4 6.5];
crystal_ls_jtou = [-53.8 0.0 0.8 18.7];
spend_ls_jtou = [-55.8 -0.8 3.5 21.0];
two_ls_jtou = [-41.6 -1.2 -0.1 7.8];
just_ls_jtou = [-44.3 0.1 1.8 18.2];

figure;
subplot(2,2,1);
y = [0 crystal_real_jtou(2:end) - crystal_ls_jtou(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 2, JTOU: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Crystal: C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,2);
y = [0 spend_real_jtou(2:end) - spend_ls_jtou(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Spend: P');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,3);
y = [0 two_real_jtou(2:end) - two_ls_jtou(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Two: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,4);
y = [0 just_real_jtou(2:end) - just_ls_jtou(2:end)];
z = [z; y(2) y(4) y(1) y(3)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Just: J');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});


% Overall
figure;
subplot(3,1,1);
y = [0 overall_real_eil(2:end) - overall_ls_eil(2:end)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('B) Subject 2: Overall Source Level Variation','fontweight','bold','fontsize',12);
xlabel('EIL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC', '', 'NG14', '', 'LG14', '', 'RR'});

subplot(3,1,2);
y = [0 overall_real_atr(2:end) - overall_ls_atr(2:end)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('ATR');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC', '', 'NG14', '', 'LG14', '', 'RR'});

subplot(3,1,3);
y = [0 overall_real_jtou(2:end) - overall_ls_jtou(2:end)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('JTOU');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'AC', '', 'NG14', '', 'LG14', '', 'RR'});


% Overall Reordered
figure;
subplot(3,1,1);
y = [overall_real_eil(2) - overall_ls_eil(2) overall_real_eil(4) - overall_ls_eil(4) 0 overall_real_eil(3) - overall_ls_eil(3)];
x = [y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('D) Subject 2: Rooms Reordered','fontweight','bold','fontsize',12);
xlabel('EIL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

subplot(3,1,2);
y = [overall_real_atr(2) - overall_ls_atr(2) overall_real_atr(4) - overall_ls_atr(4) 0 overall_real_atr(3) - overall_ls_atr(3)];
x = [x; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('ATR');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

subplot(3,1,3);
y = [overall_real_jtou(2) - overall_ls_jtou(2) overall_real_jtou(4) - overall_ls_jtou(4) 0 overall_real_jtou(3) - overall_ls_jtou(3)];
x = [x; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('JTOU');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

% Overall room only
figure;
subplot(3,1,1);
y = overall_ls_eil(2:end);
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
title('Sub. 2 - Level Relative to AC: Room Only','fontweight','bold','fontsize',12);
xlabel('EIL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});
subplot(3,1,2);
y = overall_ls_atr(2:end);
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
xlabel('ATR');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});
subplot(3,1,3);
y = overall_ls_jtou(2:end);
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
xlabel('JTOU');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});


% Averages
figure;
errorbar(mean(z), std(z), 'o','MarkerFaceColor','b');
hold on;
plot(median(z),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('B) Subject 2: Consonant Level Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Level (dB)');
hline = refline([0 0]);
set(hline,'Color','k');

figure;
errorbar(mean(x), std(x), 'o','MarkerFaceColor','b');
hold on;
plot(median(x),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('B) Subject 2: Overall Level Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Level (dB)');
hline = refline([0 0]);
set(hline,'Color','k');

% Sub 3
% EIL
overall_real_eil = [0.2 8.2 -30.4 1.2];
trees_real_eil = [1 16.7 -41.5 0.8];
hobble_real_eil = [11.6 20.6 -40.7 7.5];
crocus_real_eil = [5 11.9 -35.7 2.5];
spring_real_eil = [18.1 21.8 -58.9 15.6];

overall_ls_eil = [1 9.6 -36 1.3];
trees_ls_eil = [0.3 19 -45.4 0.9];
hobble_ls_eil = [-0.9 16.9 -43.5 2.6];
crocus_ls_eil = [-1.6 20.7 -48.5 3];
spring_ls_eil = [0.9 27.4 -63.7 6.2];

figure;
subplot(2,2,1);
y = [trees_real_eil(1) - trees_ls_eil(1) trees_real_eil(2) - trees_ls_eil(2) 0 trees_real_eil(4) - trees_ls_eil(4)];
z = [y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 3, EIL: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Trees: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,2);
y = [hobble_real_eil(1) - hobble_ls_eil(1) hobble_real_eil(2) - hobble_ls_eil(2) 0 hobble_real_eil(4) - hobble_ls_eil(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Hobble: B');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,3);
y = [crocus_real_eil(1) - crocus_ls_eil(1) crocus_real_eil(2) - crocus_ls_eil(2) 0 crocus_real_eil(4) - crocus_ls_eil(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Crocus: 2nd C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,4);
y = [spring_real_eil(1) - spring_ls_eil(1) spring_real_eil(2) - spring_ls_eil(2) 0 spring_real_eil(4) - spring_ls_eil(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Spring: P');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});


% ATR
overall_real_atr = [-1.3 8.5 -28.2 2.1];
trod_real_atr = [-4.4 23.4 -45.1 5.3];
crystal_real_atr = [-0.4 21.1 -51 3];
god_real_atr = [-5.6 17.8 -42.9 -2.9];
beautiful_real_atr = [-9.1 19.5 -36.6 6.2];

overall_ls_atr = [1.1 10.5 -34.2 1.4];
trod_ls_atr = [2 19.7 -49.2 3.1];
crystal_ls_atr = [1.6 36 -58.7 5];
god_ls_atr = [0.1 21 -48.1 0.5];
beautiful_ls_atr = [-0.8 15.6 -41.5 -0.2];

figure;
subplot(2,2,1);
y = [trod_real_atr(1) - trod_ls_atr(1) trod_real_atr(2) - trod_ls_atr(2) 0 trod_real_atr(4) - trod_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 3, ATR: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Trod: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,2);
y = [crystal_real_atr(1) - crystal_ls_atr(1) crystal_real_atr(2) - crystal_ls_atr(2) 0 crystal_real_atr(4) - crystal_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Crystal: C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,3);
y = [god_real_atr(1) - god_ls_atr(1) god_real_atr(2) - god_ls_atr(2) 0 god_real_atr(4) - god_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('God: G');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,4);
y = [beautiful_real_atr(1) - beautiful_ls_atr(1) beautiful_real_atr(2) - beautiful_ls_atr(2) 0 beautiful_real_atr(4) - beautiful_ls_atr(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Beautiful: B');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

% WDC
overall_real_wdc = [0.7 9.1 -30.9 0.4];
picture_real_wdc = [14.3 23.6 -50.9 11.5];
covers_real_wdc = [5.9 22.2 -44.8 -2.1];
darling_real_wdc = [-0.1 8 -33.4 -6.2];
violets_real_wdc = [-2.3 21.1 -42.2 -0.1];

overall_ls_wdc = [1.2 9 -36.1 0.7];
picture_ls_wdc = [0.1 26.5 -57.5 1.4];
covers_ls_wdc = [0.4 19.8 -49.1 1.1];
darling_ls_wdc = [-0.4 11.3 -40.5 0.2];
violets_ls_wdc = [0.3 17.4 -42.8 0.8];

figure;
subplot(2,2,1);
y = [picture_real_wdc(1) - picture_ls_wdc(1) picture_real_wdc(2) - picture_ls_wdc(2) 0 picture_real_wdc(4) - picture_ls_wdc(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('Subject 3, WDC: Consonant Source Level','fontweight','bold','fontsize',12);
xlabel('Picture: P');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,2);
y = [covers_real_wdc(1) - covers_ls_wdc(1) covers_real_wdc(2) - covers_ls_wdc(2) 0 covers_real_wdc(4) - covers_ls_wdc(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Covers: C');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,3);
y = [darling_real_wdc(1) - darling_ls_wdc(1) darling_real_wdc(2) - darling_ls_wdc(2) 0 darling_real_wdc(4) - darling_ls_wdc(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Darling: D');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

subplot(2,2,4);
y = [violets_real_wdc(1) - violets_ls_wdc(1) violets_real_wdc(2) - violets_ls_wdc(2) 0 violets_real_wdc(4) - violets_ls_wdc(4)];
z = [z; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('Violets: T');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14','RR','AC', 'LG14'});

% Overall
figure;
subplot(3,1,1);
y = [overall_real_eil(1) - overall_ls_eil(1) overall_real_eil(2) - overall_ls_eil(2) 0 overall_real_eil(4) - overall_ls_eil(4)];
x = [y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
title('C) Subject 3: Overall Source Level Variation','fontweight','bold','fontsize',12);
xlabel('EIL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

subplot(3,1,2);
y = [overall_real_atr(1) - overall_ls_atr(1) overall_real_atr(2) - overall_ls_atr(2) 0 overall_real_atr(4) - overall_ls_atr(4)];
x = [x; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('ATR');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});

subplot(3,1,3);
y = [overall_real_wdc(1) - overall_ls_wdc(1) overall_real_wdc(2) - overall_ls_wdc(2) 0 overall_real_wdc(4) - overall_ls_wdc(4)];
x = [x; y];
plot(y,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
hline = refline([0 0]);
set(hline,'Color','k');
xlabel('WDC');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', 'RR', '', 'AC', '', 'LG14'});


% Overall room only
figure;
subplot(3,1,1);
y = [overall_ls_eil(1) overall_ls_eil(4) overall_ls_eil(2)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
title('Sub. 3 - Level Relative to AC: Room Only','fontweight','bold','fontsize',12);
xlabel('EIL');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});
subplot(3,1,2);
y = [overall_ls_atr(1) overall_ls_atr(4) overall_ls_atr(2)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
xlabel('ATR');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});
subplot(3,1,3);
y = [overall_ls_wdc(1) overall_ls_wdc(4) overall_ls_wdc(2)];
plot(y,'o','MarkerFaceColor','b');
for i = 1:3
    text(i+0.1,y(i),num2str(y(i)),'Color','red','FontSize',12);
end
xlabel('WDC');
ylabel('dB change vs. AC');
set(gca,'XTickLabel',{'NG14', '', '', '', '', 'LG14', '', '', '', '', 'RR'});


% Averages
figure;
errorbar(mean(z), std(z), 'o','MarkerFaceColor','b');
hold on;
plot(median(z),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('C) Subject 3: Consonant Level Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Level (dB)');
hline = refline([0 0]);
set(hline,'Color','k');

figure;
errorbar(mean(x), std(x), 'o','MarkerFaceColor','b');
hold on;
plot(median(x),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('C) Subject 3: Overall Level Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Level (dB)');
hline = refline([0 0]);
set(hline,'Color','k');