function buildnall
load pos

for i=1:1
    i
    fname = strcat('fall',num2str(i));
    load(fname); 
    fallp = fall;
    fallr = fall;

    fname = strcat('nall',num2str(i));
    nallp = buildplacefield2(fallp,p);
    nallr = buildplacefield2(fallr,p);
    save(fname,'nallp','fallp'); 
end
