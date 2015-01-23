function out = vlc_channel_filter(in, channel)
out = filter(channel,1,in);
end
