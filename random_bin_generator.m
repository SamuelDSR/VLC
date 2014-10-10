function out = random_bin_generator(nBits)
%      s = RandStream('mt19937ar','Seed','shuffle');
%      RandStream.setGlobalStream(s);
out = randi([0 1],nBits,1);
end