function [Gs, log_likes, trials, flips] = MetroGibbs_binary(f,  L, varargin)

d = f.d;

% K is the number of binary flip proposals between recorded d-dimensional samples 

if length(varargin) >0
    K = varargin{1};
else
    K = d;
end


%S = sign(rand(d,1)-.5);
S= ones(d,1);


Gs = zeros(d,L);
Gs(:,1)=S;    

log_likes= zeros(L,1);
log_likes(1) = f.logp(S);

indices = randsample(d,L*K,true);

trials=0;
flips =0;

for i=2:L
    i
    for k = 1:K
        
        j=indices((i-2)*K+k);
        trials = trials +1;        
        qi = exp(-S(j)*f.logp_change(S,j));
        
        if rand() < min(1,qi)               
            S(j) = -S(j);    % flip it
            flips= flips +1;
        end
                        
    end
    ll = f.logp(S);
    log_likes(i) = ll;    
    Gs(:,i)=S;    
    mean(S);
     
    
end
