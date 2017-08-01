% %%%%%%%%%%%%%Growing Duplex Model (linear  Kernel) %%%%%%%%%%%%%%%%%%%
% This code generates a  undirected unweighted duplex or triplex network with 
% Poisson distribution of multilinks
%
% INPUTS: 
%
% N total number of nodes in the multiplex
% parameter 0<=a<=1
% parameter 0<=b<=1
% initial number of links m
% The output is  a cell array A of dimension 2
% A{1} is the adjacency matrix of the first layer 
%          
% A{2} is the adjacency matrix of the second layer 
%
% This code can be redistributed and/or modified
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%  
% This program is distributed ny the authors in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%  
% If you use this code please cite 
%
% [1] V.Nicosia, G. Bianconi, V. Latora, M. Barthelemy
% Growing multiplex networks
% Physical Review Letters 111, no. 5 (2013): 058701.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function [A] = Growing_Duplex_Linear(N,a,b,m)

%Initial condition
A{1}=sparse(N,N);
A{2}=sparse(N,N);
for i=1:m,
    for j=i+1:m,
        A{1}(i,j)=1;
        A{1}(j,i)=1;
    end
end
A{2}=A{1};

for i=(m+1):N,
    Z(1,:)=a*sum(A{1})+(1-a)*sum(A{2});
    Z(2,:)=(1-b)*sum(A{2})+(b)*sum(A{1});
    for n=1:2,
        occ=zeros(1,N);
        mx=0;
        nx=1;
       while(mx<m), 
            x=sum(Z(n,:))*rand(1);
             for ni=1:(i-1),
                 x=x-Z(n,ni);
                 if(x<0)
                     nx=ni;
                     break
                 end
             end
         while (occ(nx)==1)    
             x=sum(Z(n,:))*rand(1);
             for ni=1:(i-1),
                 x=x-Z(n,ni);
                 if(x<0)
                     nx=ni;
                     break
                 end
             end
         end
         mx=mx+1;
         occ(nx)=1;
         A{n}(i,nx)=1;
         A{n}(nx,i)=1;
         A{n}(i,i)=0;
       end
    end
 end
end

