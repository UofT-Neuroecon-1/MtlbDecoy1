ts = ftp('ftp.remidaviet.com','MatlabData@davserv.net','t6ubJ3Mn1qQm7gC9AAJb')
cd(ts,'Current')
mget(ts,'*.m')
close(ts);

if exist('Start_Experiment_Here.m','file')
    delete 'Start_Experiment_Here.m';
end