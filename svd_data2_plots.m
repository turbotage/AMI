%% Singular value interval starts and lengths for parametric sweep
singular_starts = [0,0.001,0.005,0.01,0.02,0.04,0.07,0.1,0.15,0.20,0.25,0.30,0.35,0.40,...
    0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95];
inters = [0,1,2,5,10,20,50,100,300,500];
%%
ratios = load('ratios_1555.mat').ratios;
ratios_tgc = load('ratios_1555_tgc.mat').ratios;

%%
createfigure(singular_starts,inters,ratios');
%%
createfigure(singular_starts,inters,ratios_tgc');