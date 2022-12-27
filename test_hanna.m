%% T-test dataset1 for red and blue channel
% medelstandardavvikelsen av standardavvikelsen per line

% rest
for i = 1:nodepths
    stdlinehir(i) = std(Bmodeshir(i,:));
    stdlinelor(i) = std(Bmodeslor(i,:));
end

%work
for i = 1:nodepths
    stdlinehi(i) = std(Bmodeshi(i,:));
    stdlinelo(i) = std(Bmodeslo(i,:));
end


% std per line at rest vs work and with highpass or lowpass
figure(8)
hold on;
h = histfit(stdlinehi); 
set(h(2),'color','g')
histfit(stdlinehir)
title('Std highpass at work vs rest')
xlabel('Std per line')
ylabel('Intesity')
legend('work', 'norm. dist.', 'rest', 'norm. dist.')

figure(9) 
hold on;
h = histfit(stdlinelo); 
set(h(2),'color','g')
histfit(stdlinelor)
title('Std lowpass at work vs rest')
xlabel('Std per line')
ylabel('Intensity')
legend('work', 'norm. dist.', 'rest', 'norm. dist.')


[h_hi,p_hi, ci_hi, stats_hi] = ttest2(stdlinehir,stdlinehi,'Alpha', 0.05, 'Vartype','unequal');
[h_lo,p_lo, ci_lo, stats_lo] = ttest2(stdlinelor,stdlinelo,'Alpha', 0.05,'Vartype','unequal');