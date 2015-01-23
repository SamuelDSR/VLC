function out = vlc_pd_filter(in, pd)
out = filter(pd,1,in);
end