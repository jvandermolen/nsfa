addpath('utils');
load example.mat
Ycentered=Y-repmat(mean(Y,2),1,size(Y,2));
Ynorm=Ycentered./repmat(std(Ycentered,0,2),1,size(Ycentered,2));
settings=defaultsettings();
[settings.D,settings.N]=size(Ynorm);
settings.iterations=100;
settings.store_samples=100;
%settings.output='ex3r3.mat';
mvmask=binornd(1,1,settings.D,settings.N);
initialsample=init_nsfa(settings);
[finalsample,resultstable]=nsfa(Ynorm,mvmask,initialsample,settings);
plot(resultstable(:,1),resultstable(:,2),'k+-'); 
xlabel('cpu time');
ylabel('log joint probability');
