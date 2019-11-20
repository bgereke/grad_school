function virmenAboutWindow

fig = figure;
backupUnits = get(0,'units');
set(0,'units','pixels');
scr = get(0,'screensize');
wd = 500;
hg = 200;
set(fig,'position',[scr(3)/2-wd/2 scr(4)/2-hg/2 wd hg],'resize','off', ...
    'menubar','none','name','About ViRMEn','numbertitle','off','windowstyle','modal');
set(0,'units',backupUnits);

txt = text(0,1,'{\color{red}{\bfVi}}rtual {\color{red}{\bfR}}eality {\color{red}{\bfM}}atlab {\color{red}{\bfEn}}gine');
set(txt,'horizontalalignment','center','verticalalignment','top','fontsize',20);

txt = text(0,.45,'by Dmitriy Aronov, 2011');
set(txt,'horizontalalignment','center','verticalalignment','middle','fontsize',12);

mfile = mfilename('fullpath');
path = fileparts(mfile);
load([path filesep 'virmenVersion.mat']);
txt = text(0,.15,['{\bfCurrent version:} ' virmenVersion]);
set(txt,'horizontalalignment','center','verticalalignment','middle','fontsize',12);

xlim([-1 1]);
ylim([0 1]);
axis off