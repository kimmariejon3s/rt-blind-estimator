clear all;
close all;
f = [63 125 250 500 1000 2000 4000 8000];
balloon_t20 = [0.753568613 0.734623896 0.77602034 0.710522167 0.580505205 0.481879427 0.473978709 0.424601904];
balloon_se = [0.012449677 0.023828171 0.016232584 0.008242777 0.006130813 0.005302714 0.002502875 0.003423653];
balloon_se = 1.96 .* balloon_se;        % 95% confidence interval

annie_3min = [
    NaN NaN;
    NaN NaN;
    0.602200787421647   0.059657683475808; 
    0.586258560127753   0.026863588783899;
    0.524333576854855   0.030709224064383;
    0.451574396114073   0.018549813064614;
    0.474617692093611   0.016280917449663;
    0.415831872905924   0.013892097493394];

%pos = ((balloon_t20 - balloon_se) + (annie_3min(:,1) + annie_3min(:,2))') ./ 2;
limen = balloon_t20 .* 0.05;

annie_6min = [
    NaN NaN;
    NaN NaN;
    0.530317670266028   0.090622206699213;
    0.556043119782238   0.038550355124205;
    0.499249654333680   0.041629292327950;
    0.438347779669901   0.023550267019663;
    0.467021092805973   0.021191495443792;
    0.401413926971080   0.013984204766845];

summertime_3min = [
    NaN   NaN;
    NaN   NaN;
    NaN   NaN;
    NaN   NaN;
    0.604054706950595   0.039290399280173;
    NaN   NaN;
    NaN   NaN;
    NaN   NaN];

summertime_6min = [
    NaN   NaN;
    NaN   NaN;
    0.772715772798639   0.077197975928000;
    NaN NaN;
    0.595445969950027   0.036062858757722;
    0.561607478598799   0.042684378744356;
    0.585898165055409   0.054842335130246;
    NaN NaN];

both_3min = [
    NaN   NaN;
    NaN   NaN;
    0.652407853922697   0.069154861554943;
    0.565313718404708   0.036734986522291;
    0.521416418657599   0.034278457014846;
    0.465240440843934   0.027952469434685;
    0.472992882572143   0.013631131990024;
    0.406859605575855   0.014443932530305;
];

both_6min = [
    NaN   NaN;
    NaN   NaN;
    0.687674922543754   0.067264309250359;
    0.603325259733473   0.072100732267561;
    0.504938456881736   0.037023829614642;
    0.483208769504158   0.037373552678454;
    0.541266126896288   0.061309042974778;
    NaN NaN;
];

errorbar(balloon_t20, balloon_se,'b');
hold on;
errorbar(annie_3min(:,1), annie_3min(:,2),'g');
errorbar(summertime_3min(:,1), summertime_3min(:,2),'r');
%errorbar(both_3min(:,1), both_3min(:,2),'k');

title('RT60 (3 min segments)');
legend('Balloon', 'Hard Knock', 'Summertime','Both Tracks','Location','Northeast');
ylabel('T20 (s)');
xlabel('Octave Band Centre Frequency (Hz)');
set(gca,'XTickLabel',{'0','63','125','250','500','1000','2000','4000','8000','16000'});

figure(2);
errorbar(balloon_t20, balloon_se,'b');
hold on;
%errorbar(annie_6min(:,1), annie_6min(:,2),'g');
errorbar(summertime_3min(:,1), summertime_3min(:,2),'g');
errorbar(summertime_6min(:,1), summertime_6min(:,2),'r');
%errorbar(both_3min(:,1), both_3min(:,2),'k');

title('RT60 (Summertime)');
legend('Balloon', 'Summertime 3 min', 'Summertime 6 min','Location','Northeast');
ylabel('T20 (s)');
xlabel('Octave Band Centre Frequency (Hz)');
set(gca,'XTickLabel',{'0','63','125','250','500','1000','2000','4000','8000','16000'});

figure(3);
errorbar(annie_3min(:,1), annie_3min(:,2),'r');
hold on;
errorbar(annie_6min(:,1), annie_6min(:,2),'g');
title('RT60 (Hard Knock Life)');
legend('3 min windows', '6 min windows', 'Location','Northeast');
ylabel('T20 (s)');
xlabel('Octave Band Centre Frequency (Hz)');
set(gca,'XTickLabel',{'0','63','125','250','500','1000','2000','4000','8000','16000'});

figure(4);
errorbar(summertime_3min(:,1), summertime_3min(:,2),'r');
hold on;
errorbar(summertime_6min(:,1), summertime_6min(:,2),'g');
title('RT60 (Summertime)');
legend('3 min windows', '6 min windows', 'Location','Northeast');
ylabel('T20 (s)');
xlabel('Octave Band Centre Frequency (Hz)');
set(gca,'XTickLabel',{'0','63','125','250','500','1000','2000','4000','8000','16000'});

figure;
errorbar(balloon_t20, limen,'b');
hold on;
errorbar(annie_3min(:,1), annie_3min(:,2),'g');
title('RT60 of Staccato Singing (3 min segments) showing Difference Limen (DL)');
legend('Balloon w. DL', 'Hard Knock','Location','Northeast');
ylabel('T20 (s)');
xlabel('Octave Band Centre Frequency (Hz)');
set(gca,'XTickLabel',{'0','63','125','250','500','1000','2000','4000','8000','16000'});