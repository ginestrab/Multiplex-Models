% %%%%%%%%%%%%%Growing Duplex Model (Non-linear  Kernel) %%%%%%%%%%%%%%%%%%%
% This code generates a  undirected unweighted duplex or triplex network with 
% Poisson distribution of multilinks
%
% INPUTS: 
%
% N total number of nodes in the multiplex
% parameter alpha 
% parameter beta
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
%  "Nonlinear growth and condensation in multiplex networks." 
% Physical Review E 90, no. 4 (2014): 042807.0
% 
function [A] = Growing_Duplex_Nonlinear(N,alpha,beta,m)


%Initial condition
A{1}=sparse(N,N);
A{2}=sparse(N,N);
for i=1:m;
    for j=i+1:m;
        A{1}(i,j)=1;
        A{1}(j,i)=1;
    end
end
A{1}(m,m)=0;
A{2}=A{1};

for i=(m+1):N,
    

    y1=sum(A{1});
    y2=sum(A{2});
    if(abs(alpha*beta)>0)
    Z(1,:)=y1.*(y1+(y1==0)).^(alpha-1).*y2.*(y2+(y2==0)).^(beta-1);
    Z(2,:)=y2.*(y2+(y2==0)).^(alpha-1).*y1.*(y1+(y1==0)).^(beta-1);
    end
    if((alpha==0)&&(abs(beta)>0)),
    Z(1,:)=y2.*(y2+(y2==0)).^(beta-1);
    Z(2,:)=y1.*(y1+(y1==0)).^(beta-1);
    end
    if((beta==0)&&(abs(alpha)>0)),
        Z(1,:)=y1.*(y1+(y1==0)).^(alpha-1);
        Z(2,:)=y2.*(y2+(y2==0)).^(alpha-1);
    end
    if((beta==0)&&(alpha==0)),
        Z(1,:)=(y2+(y2==0)).^(alpha);
        Z(2,:)=(y2+(y2==0)).^(alpha);
    end
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
                     break;
                 end
             end
         while (occ(nx)==1)    
             x=sum(Z(n,:))*rand(1);
             for ni=1:(i-1),
                 x=x-Z(n,ni);
                 if(x<0)
                     nx=ni;
                     break;
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