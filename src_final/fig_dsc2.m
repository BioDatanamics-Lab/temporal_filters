clear all; tic; close all;

axFnt = 12; ttlFnt = 23; lblFnt = 14;

ODEparams.E_L = -60; ODEparams.E_ex = 0; ODEparams.C = 1; ODEparams.U = .1;
ODEparams.G_ex = 1;  ODEparams.I_app = 0; ODEparams.tau_dec = 3;

f1 = figure(1); set(gcf,'Position', [10 10 900 500]);

gL = [.08 .1]; dep = [20 400];
ODEparams.tau_fac = 0;
reps = 1100;
for kkk = 1:2
for kk = 1:2
ODEparams.g_L = gL(kk);
ODEparams.tau_dep = dep(kkk);

%%% regiemes will be a N by 2 matrix, rows will be different reigiems,
%%% columns will be: 1 - Hz, 2- ms of simulation for that window
regiemes = [20 7500; 40 5000; 80 2000];
dmy = sum(regiemes,1);
storage = zeros(1,dmy(2)/.01 + 1);
[~, timecourse, ~] = numerical_dscPois(ODEparams,regiemes);

parfor i = 1:reps
   tic
   [soln, timecourse, isFired] = numerical_dscPois(ODEparams,regiemes);
   storage = squeeze(soln(4,:))+storage;
   toc
end

if kkk == 1 && kk == 1
plot(timecourse,storage/reps,'-b'); hold on;
1
elseif kkk == 1 && kk == 2
    plot(timecourse,storage/reps,'Color',[0.5843 0.8157 0.9882]); hold on;
    2
elseif kkk == 2 && kk == 1
plot(timecourse,storage/reps,'-r'); hold on;
3
elseif kkk == 2 && kk == 2
    plot(timecourse,storage/reps,'Color',[255 153 153]/256); hold on;
    4
end

end
end

xlim([3000 14500])
legend('G_L = .08, \tau_{dep} = 10','G_L = 1, \tau_{dep} = 10','G_L = .08, \tau_{dep} = 250','G_L = 1, \tau_{dep} = 250','Location','northwest')
ylabel('V [mV]'); xlabel('t [ms]')

%saveas(f1,'../../img/TF/fig_dsc2','epsc')
