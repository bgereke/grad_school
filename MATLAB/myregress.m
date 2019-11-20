function out = myregress(x,y)

X = [ones(size(x,2),1) x(1,:)'];
y = y';
[B,~,~,~,STATS] = regress(y(:,randsample(size(y,1),1)),X);
out = [B(2),STATS];
end