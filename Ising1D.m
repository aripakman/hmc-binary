 classdef Ising1D < handle  % inherit from handle so that we can pass by reference

        % this class implements a 1d Ising model with zero magnetic field. 
     
        properties(SetObservable = true)
        % These properties are public by default
           d;
           beta;
           mlp=-Inf;   
           M;
           Neis;
        end
        
             
        methods
        % These methods are public by default. 
        
            function obj = Ising1D(d,T)
            % class constructor            
            
                obj.d = d;    %linear dimension of the 2D grid
                obj.beta=1/T;
                
                obj.M = sparse(obj.d,obj.d);
                obj.Neis = zeros(obj.d,2);                
                
                    for j=1:d                            
                        nei = neighbors(j);
                        obj.Neis(j,:) = nei;                            
                        obj.M(j, nei) =-1;                
                    end                                
                    obj.M = obj.beta*obj.M/2;    % the factor of 1/2 is because pairs are counted twice                

                function nei = neighbors(j)
                    if j == 1
                        nei = [d, 2];
                    elseif j == d
                        nei = [d-1,1  ]; 
                    else
                        nei = [j-1, j+1];
                    end                        
                end                
            end
            
            
                
            
            % Conventions:
            % S is a signs vector  (+1,-1)
            % H = -log P(S) = S'*obj.M*S 
            
             function lp = logp(obj,S)                             
                lp = -S'*obj.M*S;
                if lp > obj.mlp
                    obj.mlp = lp;
                end                
                    
             end
         
             function lpc= logp_change(obj,S,j)
                 % returns the difference in the log probability when 
                 % S(j) == +1 and S(j) == -1     
                 
                 nei = obj.Neis(j,:);
                 lpc = 2*obj.beta*sum(S(nei));                 
                     
             end
             
             
 
             
        end
        
        
        
end
