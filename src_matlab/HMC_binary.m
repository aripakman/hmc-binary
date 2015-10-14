function [Zs, log_likes, energies] = HMC_binary(f, T,L)


% Implementation of the Gaussian augmented-variable HMC sampler introduced in the NIPS 2013 paper 
% "Auxiliary-variable exact Hamiltonian Monte Carlo samplers for binary distributions" by Ari Pakman and Liam Paninski 
% Author: Ari Pakman

% Input:
% f:          object that contains the binary distribution of interest. 
%             It implements the methods f.d, f.logp_change(S,j) and f.logp(S)
% L:          number of samples desired
% T:          travel time for the HMC

% Output
% Zs:      d x L matrix, each column is a sample

%% 


d = f.d;

initial_X = abs(normrnd(0,1,d,1));

log_likes = zeros(L,1);
energies  = zeros(L,1);

ll = f.logp(sign(initial_X));
log_likes(1) = ll;        


Xs=NaN;

wall_hits =0;
wall_crosses =0;


nearzero= 10000*eps;

touched = zeros(d,int64(T/pi)*(L-1));
mts = zeros(d,int64(T/pi)*(L-1));
xx = zeros(d,1);


%% Sampling loop


last_X= initial_X;
Xs=zeros(d,L);


Xs(:,1)=initial_X;


i=2;
while (i <= L)
i
stop=0;   
j=0;
V= normrnd(0,1, d,1);   % initial velocity
X = last_X;

tt=0;                    % records how much time the particle already moved 
S=sign(X);



    while (1)
        
        a = V; 
        b = X;        
        phi = atan2(b,a);           % -pi < phi < +pi    

        % find the first time constraint becomes zero

            wt1= -phi;                 % time at which coordinates hit the walls                                             % wt1 in [-pi/2, 3pi/3]
            wt1(phi>0) = pi -phi(phi>0); 


         % if there was a previous reflection (j>0)
         % and there is a potential reflection at the sample plane                                    
         % make sure that a new reflection at j is not found because of numerical error
         
            if j>0    
                    tt1 = wt1(j);
                    if abs(tt1) < nearzero || abs(tt1-2*pi)< nearzero
                        wt1(j)=Inf;
                    end                    
            end


            [mt, j] = min(wt1);

        tt=tt+mt;
%       fprintf(num2str(tt/T));

        if tt>=T
            mt= mt-(tt-T);
            stop=1;

        else
            wall_hits = wall_hits + 1;
            

        end

        % move the particle a time mt

        X = a*sin(mt) + b*cos(mt);
        V = a*cos(mt) - b*sin(mt);

        if stop                    
            break;
        end

        X(j) = 0;

        v2_new = V(j)^2 +  sign(V(j))*2*f.logp_change(S,j);
        if v2_new >0
            V(j) = sqrt(v2_new)* sign(V(j));
            S(j) = -S(j);
            wall_crosses = wall_crosses +1;
   %        disp('crossed');

        else
            V(j) = -V(j);
        end
        
        
        
    end % while(1)

    % at this point we have a sampled value X
    
        Xs(:,i)=X;
        ll = f.logp(S);
        log_likes(i) = ll;                
        energies(i) = 0.5*(X'*X + V'*V) -log_likes(i);        
        last_X = X;
        i= i+1;
        mean(S);

        
end %while (i <= L)


Sis= sign(Xs);

% to transform the sampled values into binary vectors, with 0's and 1's.
% Zs= (Sis + ones(size(Sis)))/2;
Zs = Sis;


end

