clear all;
close all;

%%% tempo sub 1
sub1_eil = [84 81 86 81];
sub1_atr = [57 53 54 54];
sub1_wndl = [79 71 78 79];

figure;
subplot(3,1,1);
plot(sub1_eil,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub1_eil(i),num2str(sub1_eil(i)),'Color','red','FontSize',12);
end
title('A)   Subject 1: Tempo','fontweight','bold','fontsize',12);
ylabel('BPM');
xlabel('Everywhere I Look');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

subplot(3,1,2);
plot(sub1_atr,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub1_atr(i),num2str(sub1_atr(i)),'Color','red','FontSize',12);
end
xlabel('At The River');
ylabel('BPM');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

subplot(3,1,3);
plot(sub1_wndl,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub1_wndl(i),num2str(sub1_wndl(i)),'Color','red','FontSize',12);
end
xlabel('Wer Nur Den Lieben Gott Lasst Walten');
ylabel('BPM');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

%%% tempo sub 2
sub2_eil = [107 107 103 100];
sub2_atr = [58 64 64 60];
sub2_jtou = [98 101 104 103.5];

figure;
subplot(3,1,1);
plot(sub2_eil,'o','MarkerFaceColor','b')
for i = 1:4
    text(i+0.1,sub2_eil(i),num2str(sub2_eil(i)),'Color','red','FontSize',12);
end
title('B)   Subject 2: Tempo','fontweight','bold','fontsize',12);
ylabel('BPM');
xlabel('Everywhere I Look');
set(gca,'XTickLabel',{'AC','','NG14','','LG14','','RR'});

subplot(3,1,2);
plot(sub2_atr,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub2_atr(i),num2str(sub2_atr(i)),'Color','red','FontSize',12);
end
xlabel('At The River');
ylabel('BPM');
set(gca,'XTickLabel',{'AC','','NG14','','LG14','','RR'});

subplot(3,1,3);
plot(sub2_jtou,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub2_jtou(i),num2str(sub2_jtou(i)),'Color','red','FontSize',12);
end
xlabel('Just the Two of Us');
ylabel('BPM');
set(gca,'XTickLabel',{'AC','','NG14','','LG14','','RR'});

%%% tempo sub 3
sub3_eil = [111.5 108.5 110 108.5];
sub3_atr = [60 59.5 59.5 59.5];
sub3_wdc = [128 126 128 128];

figure;
subplot(3,1,1);
plot(sub3_eil,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub3_eil(i),num2str(sub3_eil(i)),'Color','red','FontSize',12);
end
title('C)   Subject 3: Tempo','fontweight','bold','fontsize',12);
ylabel('BPM');
xlabel('Everywhere I Look');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

subplot(3,1,2);
plot(sub3_atr,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub3_atr(i),num2str(sub3_atr(i)),'Color','red','FontSize',12);
end
xlabel('At The River');
ylabel('BPM');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

subplot(3,1,3);
plot(sub3_wdc,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub3_wdc(i),num2str(sub3_wdc(i)),'Color','red','FontSize',12);
end
xlabel('When Doves Cry');
ylabel('BPM');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});


%%% tempo sub 2 same order as 1 and 3
sub2_eil = [107 100 107 103];
sub2_atr = [64 60 58 64];
sub2_jtou = [101 103.5 98 104];

figure;
subplot(3,1,1);
plot(sub2_eil,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub2_eil(i),num2str(sub2_eil(i)),'Color','red','FontSize',12);
end
title('D)   Subject 2: Tempo (rooms reordered)','fontweight','bold','fontsize',12);
ylabel('BPM');
xlabel('Everywhere I Look');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

subplot(3,1,2);
plot(sub2_atr,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub2_atr(i),num2str(sub2_atr(i)),'Color','red','FontSize',12);
end
xlabel('At The River');
ylabel('BPM');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

subplot(3,1,3);
plot(sub2_jtou,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,sub2_jtou(i),num2str(sub2_jtou(i)),'Color','red','FontSize',12);
end
xlabel('Just the Two of Us');
ylabel('BPM');
set(gca,'XTickLabel',{'NG14','','RR','','AC','','LG14'});

% Average relative tempos
y = [sub1_eil - sub1_eil(3); sub1_atr - sub1_atr(3); sub1_wndl - sub1_wndl(3)];
suball_avg = mean(y);
suball_std = std(y);
figure;
errorbar(suball_avg,suball_std,'o','MarkerFaceColor','b');
hold on;
plot(median(y),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('A) Subject 1: Average Tempo Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('BPM');
hline = refline([0 0]);
set(hline,'Color','k');

y = [sub2_eil - sub2_eil(3); sub2_atr - sub2_atr(3); sub2_jtou - sub2_jtou(3)];
suball_avg = mean(y);
suball_std = std(y);
figure;
errorbar(suball_avg,suball_std,'o','MarkerFaceColor','b');
hold on;
plot(median(y),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('B) Subject 2: Average Tempo Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('BPM');
hline = refline([0 0]);
set(hline,'Color','k');

y = [ sub3_eil - sub3_eil(3); sub3_atr - sub3_atr(3); sub3_wdc - sub3_wdc(3)];
suball_avg = mean(y);
suball_std = std(y);
figure;
errorbar(suball_avg,suball_std,'o','MarkerFaceColor','b');
hold on;
plot(median(y),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('C) Subject 3: Average Tempo Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('BPM');
hline = refline([0 0]);
set(hline,'Color','k');
%%

%%% Consonant Duration Sub 1
% EIL
con_trees = [0.101 0.121 0.106 0.075];
con_hobble = [0.071 0.052 0.052 0.058];
con_crocus = [0.092 0.099 0.058 0.073];
con_spring = [0.073 0.078 0.080 0.068];
y = [con_trees - con_trees(3); con_hobble - con_hobble(3); con_crocus - con_crocus(3); con_spring - con_spring(3)];

figure;
subplot(2,2,1);
plot(con_trees,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trees(i),num2str(con_trees(i)),'Color','red','FontSize',12);
end
title('Sub. 1 Consonant Duration: EIL','fontweight','bold','fontsize',12);
xlabel('Trees: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_hobble,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_hobble(i),num2str(con_hobble(i)),'Color','red','FontSize',12);
end
xlabel('Hobble: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_crocus,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crocus(i),num2str(con_crocus(i)),'Color','red','FontSize',12);
end
xlabel('Crocus: 2nd C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_spring,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_spring(i),num2str(con_spring(i)),'Color','red','FontSize',12);
end
xlabel('Spring: P');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

%ATR
con_trod = [0.044 0.066 0.064 0.057];
con_crystal = [0.052 0.097 0.043 0.058];
con_god = [0.053 0.053 0.058 0.075];
con_beautiful = [0.084 0.078 0.078 0.089];
y = [y; con_trod - con_trod(3); con_crystal - con_crystal(3); con_god - con_god(3); con_beautiful - con_beautiful(3)];

figure;
subplot(2,2,1);
plot(con_trod,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trod(i),num2str(con_trod(i)),'Color','red','FontSize',12);
end
title('Sub. 1 Consonant Duration: ATR','fontweight','bold','fontsize',12);
xlabel('Trod: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_crystal,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crystal(i),num2str(con_crystal(i)),'Color','red','FontSize',12);
end
xlabel('Crystal: C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_god,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_god(i),num2str(con_god(i)),'Color','red','FontSize',12);
end
xlabel('God: G');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_beautiful,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_beautiful(i),num2str(con_beautiful(i)),'Color','red','FontSize',12);
end
xlabel('Beautiful: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});


%WNDL
con_den = [0.024 0.023 0.023 0.025];
con_lieben = [0.068 0.078 0.060 0.067];
con_hoffet = [0.029 0.051 0.054 0.043];
con_keit = [0.080 0.126 0.052 0.110];
y = [y; con_den - con_den(3); con_lieben - con_lieben(3); con_hoffet - con_hoffet(3); con_keit - con_keit(3)];

figure;
subplot(2,2,1);
plot(con_den,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_den(i),num2str(con_den(i)),'Color','red','FontSize',12);
end
title('Sub. 1 Consonant Duration: WNDL','fontweight','bold','fontsize',12);
xlabel('Den: D');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_lieben,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_lieben(i),num2str(con_lieben(i)),'Color','red','FontSize',12);
end
xlabel('Lieben: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_hoffet,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_hoffet(i),num2str(con_hoffet(i)),'Color','red','FontSize',12);
end
xlabel('Hoffet: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_keit,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_keit(i),num2str(con_keit(i)),'Color','red','FontSize',12);
end
xlabel('Traurigkeit: K');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

% Averages
figure;
errorbar(mean(y), std(y),'o','MarkerFaceColor','b');
hold on;
plot(median(y),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('A) Subject 1: Consonant Duration Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Duration (s)');
hline = refline([0 0]);
set(hline,'Color','k');

%%% Consonant Duration Sub 2
% EIL
con_trees = [0.047 0.047 0.057 0.044];
con_hobble = [0.038 0.032 0.052 0.070];
con_crocus = [0.058 0.041 0.045 0.042];
con_spring = [0.062 0.071 0.078 0.048];

figure;
subplot(2,2,1);
plot(con_trees,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trees(i),num2str(con_trees(i)),'Color','red','FontSize',12);
end
title('Sub. 2 Consonant Duration: EIL','fontweight','bold','fontsize',12);
xlabel('Trees: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,2);
plot(con_hobble,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_hobble(i),num2str(con_hobble(i)),'Color','red','FontSize',12);
end
xlabel('Hobble: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,3);
plot(con_crocus,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crocus(i),num2str(con_crocus(i)),'Color','red','FontSize',12);
end
xlabel('Crocus: 2nd C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,4);
plot(con_spring,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_spring(i),num2str(con_spring(i)),'Color','red','FontSize',12);
end
xlabel('Spring: P');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});


%ATR
con_trod = [0.053 0.060 0.052 0.048];
con_crystal = [0.041 0.050 0.050 0.065];
con_god = [0.061 0.060 0.053 0.071];
con_beautiful = [0.058 0.058 0.063 0.065];

figure;
subplot(2,2,1);
plot(con_trod,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trod(i),num2str(con_trod(i)),'Color','red','FontSize',12);
end
title('Sub. 2 Consonant Duration: ATR','fontweight','bold','fontsize',12);
xlabel('Trod: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,2);
plot(con_crystal,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crystal(i),num2str(con_crystal(i)),'Color','red','FontSize',12);
end

xlabel('Crystal: C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,3);
plot(con_god,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_god(i),num2str(con_god(i)),'Color','red','FontSize',12);
end
xlabel('God: G');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,4);
plot(con_beautiful,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_beautiful(i),num2str(con_beautiful(i)),'Color','red','FontSize',12);
end
xlabel('Beautiful: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});


%JTOU
con_crystal = [0.056 0.049 0.037 0.062];
con_spend = [0.080 0.077 0.085 0.107];
con_two = [0.077 0.084 0.114 0.099];
con_just = [0.070 0.075 0.073 0.082];

figure;
subplot(2,2,1);
plot(con_crystal,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crystal(i),num2str(con_crystal(i)),'Color','red','FontSize',12);
end
title('Sub. 2 Consonant Duration: JTOU','fontweight','bold','fontsize',12);
xlabel('Crystal: C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,2);
plot(con_spend,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_spend(i),num2str(con_spend(i)),'Color','red','FontSize',12);
end
xlabel('Spend: P');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,3);
plot(con_two,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_two(i),num2str(con_two(i)),'Color','red','FontSize',12);
end
xlabel('Two: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});

subplot(2,2,4);
plot(con_just,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_just(i),num2str(con_just(i)),'Color','red','FontSize',12);
end
xlabel('Just: J');
ylabel('Time (s)');
set(gca,'XTickLabel',{'AC','NG14','LG14','RR'});


%%% Consonant Duration Sub 3
% EIL
con_trees = [0.074 0.056 0.066 0.082];
con_hobble = [0.042 0.051 0.044 0.034];
con_crocus = [0.057 0.035 0.085 0.087];
con_spring = [0.082 0.091 0.065 0.053];
y = [con_trees - con_trees(3); con_hobble - con_hobble(3); con_crocus - con_crocus(3); con_spring - con_spring(3)];

figure;
subplot(2,2,1);
plot(con_trees,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trees(i),num2str(con_trees(i)),'Color','red','FontSize',12);
end
title('Sub. 3 Consonant Duration: EIL','fontweight','bold','fontsize',12);
xlabel('Trees: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_hobble,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_hobble(i),num2str(con_hobble(i)),'Color','red','FontSize',12);
end
xlabel('Hobble: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_crocus,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crocus(i),num2str(con_crocus(i)),'Color','red','FontSize',12);
end
xlabel('Crocus: 2nd C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_spring,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_spring(i),num2str(con_spring(i)),'Color','red','FontSize',12);
end
xlabel('Spring: P');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});


%ATR
con_trod = [0.063 0.058 0.034 0.030];
con_crystal = [0.055 0.082 0.057 0.058];
con_god = [0.071 0.057 0.058 0.068];
con_beautiful = [0.072 0.080 0.062 0.076];
y = [y; con_trod - con_trod(3); con_crystal - con_crystal(3); con_god - con_god(3); con_beautiful - con_beautiful(3)];

figure;
subplot(2,2,1);
plot(con_trod,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trod(i),num2str(con_trod(i)),'Color','red','FontSize',12);
end
title('Sub. 3 Consonant Duration: ATR','fontweight','bold','fontsize',12);
xlabel('Trod: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_crystal,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crystal(i),num2str(con_crystal(i)),'Color','red','FontSize',12);
end
xlabel('Crystal: C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_god,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_god(i),num2str(con_god(i)),'Color','red','FontSize',12);
end
xlabel('God: G');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_beautiful,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_beautiful(i),num2str(con_beautiful(i)),'Color','red','FontSize',12);
end
xlabel('Beautiful: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});


%WDC
con_picture = [0.185 0.161 0.147 0.124];
con_covers = [0.107 0.100 0.095 0.104];
con_darling = [0.109 0.107 0.105 0.122];
con_violets = [0.020 0.034 0.035 0.036];
y = [y; con_picture - con_picture(3); con_covers - con_covers(3); con_darling - con_darling(3); con_violets - con_violets(3)];

figure;
subplot(2,2,1);
plot(con_picture,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_picture(i),num2str(con_picture(i)),'Color','red','FontSize',12);
end
title('Sub. 3 Consonant Duration: WDC','fontweight','bold','fontsize',12);
xlabel('Picture: P');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_covers,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_covers(i),num2str(con_covers(i)),'Color','red','FontSize',12);
end
xlabel('Covers: C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_darling,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_darling(i),num2str(con_darling(i)),'Color','red','FontSize',12);
end
xlabel('Darling: D');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_violets,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_violets(i),num2str(con_violets(i)),'Color','red','FontSize',12);
end
xlabel('Violets: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

% Averages
figure;
errorbar(mean(y), std(y),'o','MarkerFaceColor','b');
hold on;
plot(median(y),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('C) Subject 3: Consonant Duration Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Duration (s)');
hline = refline([0 0]);
set(hline,'Color','k');

%%%

%%% Consonant Duration Sub 2 (rooms reordered to match sub 1 and 3)
% EIL
con_trees = [0.047 0.047 0.057 0.044];
con_hobble = [0.038 0.032 0.052 0.070];
con_crocus = [0.058 0.041 0.045 0.042];
con_spring = [0.062 0.071 0.078 0.048];

con_trees = [con_trees(2) con_trees(4) con_trees(1) con_trees(3)];
con_hobble = [con_hobble(2) con_hobble(4) con_hobble(1) con_hobble(3)];
con_crocus = [con_crocus(2) con_crocus(4) con_crocus(1) con_crocus(3)];
con_spring = [con_spring(2) con_spring(4) con_spring(1) con_spring(3)];
y = [con_trees - con_trees(3); con_hobble - con_hobble(3); con_crocus - con_crocus(3); con_spring - con_spring(3)];

figure;
subplot(2,2,1);
plot(con_trees,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trees(i),num2str(con_trees(i)),'Color','red','FontSize',12);
end
title('Sub. 2: EIL (rooms reordered)','fontweight','bold','fontsize',12);
xlabel('Trees: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_hobble,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_hobble(i),num2str(con_hobble(i)),'Color','red','FontSize',12);
end
xlabel('Hobble: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_crocus,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crocus(i),num2str(con_crocus(i)),'Color','red','FontSize',12);
end
xlabel('Crocus: 2nd C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_spring,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_spring(i),num2str(con_spring(i)),'Color','red','FontSize',12);
end
xlabel('Spring: P');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});


%ATR
con_trod = [0.053 0.060 0.052 0.048];
con_crystal = [0.041 0.050 0.050 0.065];
con_god = [0.061 0.060 0.053 0.071];
con_beautiful = [0.058 0.058 0.063 0.065];

con_trod = [con_trod(2) con_trod(4) con_trod(1) con_trod(3)];
con_crystal = [con_crystal(2) con_crystal(4) con_crystal(1) con_crystal(3)];
con_god = [con_god(2) con_god(4) con_god(1) con_god(3)];
con_beautiful = [con_beautiful(2) con_beautiful(4) con_beautiful(1) con_beautiful(3)];
y = [y; con_trod - con_trod(3); con_crystal - con_crystal(3); con_god - con_god(3); con_beautiful - con_beautiful(3)];

figure;
subplot(2,2,1);
plot(con_trod,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_trod(i),num2str(con_trod(i)),'Color','red','FontSize',12);
end
title('Sub. 2: ATR (rooms reordered)','fontweight','bold','fontsize',12);
xlabel('Trod: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_crystal,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crystal(i),num2str(con_crystal(i)),'Color','red','FontSize',12);
end
xlabel('Crystal: C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_god,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_god(i),num2str(con_god(i)),'Color','red','FontSize',12);
end
xlabel('God: G');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_beautiful,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_beautiful(i),num2str(con_beautiful(i)),'Color','red','FontSize',12);
end
xlabel('Beautiful: B');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});


%JTOU
con_crystal = [0.056 0.049 0.037 0.062];
con_spend = [0.080 0.077 0.085 0.107];
con_two = [0.077 0.084 0.114 0.099];
con_just = [0.070 0.075 0.073 0.082];

con_crystal = [con_crystal(2) con_crystal(4) con_crystal(1) con_crystal(3)];
con_spend = [con_spend(2) con_spend(4) con_spend(1) con_spend(3)];
con_two = [con_two(2) con_two(4) con_two(1) con_two(3)];
con_just = [con_just(2) con_just(4) con_just(1) con_just(3)]; 
y = [y; con_crystal - con_crystal(3); con_spend - con_spend(3); con_two - con_two(3); con_just - con_just(3)];

figure;
subplot(2,2,1);
plot(con_crystal,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_crystal(i),num2str(con_crystal(i)),'Color','red','FontSize',12);
end
title('Sub. 2: JTOU (rooms reordered)','fontweight','bold','fontsize',12);
xlabel('Crystal: C');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,2);
plot(con_spend,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_spend(i),num2str(con_spend(i)),'Color','red','FontSize',12);
end
xlabel('Spend: P');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,3);
plot(con_two,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_two(i),num2str(con_two(i)),'Color','red','FontSize',12);
end
xlabel('Two: T');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

subplot(2,2,4);
plot(con_just,'o','MarkerFaceColor','b');
for i = 1:4
    text(i+0.1,con_just(i),num2str(con_just(i)),'Color','red','FontSize',12);
end
xlabel('Just: J');
ylabel('Time (s)');
set(gca,'XTickLabel',{'NG14','RR','AC','LG14'});

% Averages
figure;
errorbar(mean(y), std(y),'o','MarkerFaceColor','b');
hold on;
plot(median(y),'ro','MarkerFaceColor','r');
legend('Mean', 'Median');
set(gca,'XTickLabel',{'', 'NG14','','RR','','AC','','LG14'});
title('B) Subject 2: Consonant Duration Differences Relative to AC','fontweight','bold','fontsize',12);
ylabel('Duration (s)');
hline = refline([0 0]);
set(hline,'Color','k');
