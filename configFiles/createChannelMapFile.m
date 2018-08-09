
load('H:\S_batch\180215S\180215\180215bin_002\mapping.mat')

Nchannels = 1024;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   =  getfield(mapping, 'x');
ycoords   = getfield(mapping, 'y');
ID=find(xcoords==-1)
connected(ID) = 0
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
fs = 20000; % sampling frequency
save('H:\S_batch\180215S\180215\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%

load('H:\S_batch\180215S\180215\180215bin_002\mapping.mat')

Nchannels = 732;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoordst   =  getfield(mapping, 'x');
ycoordst   = getfield(mapping, 'y');
xcoords=xcoordst(xcoordst~=-1)
ycoords=ycoordst(ycoordst~=-1)
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 20000; % sampling frequency
save('H:\S_batch\180215S\180215\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')