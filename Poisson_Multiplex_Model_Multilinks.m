% %%%%%%%%%%%%%Poisson Multiplex Model Multilinks %%%%%%%%%%%%%%%%%%%
% This code generates a  undirected unweighted duplex or triplex network with 
% Poisson distribution of multilinks
%
% INPUTS: 
%
% N total number of nodes in the multiplex
% M total number of layers (M=2 or M=3)
% c average degree of multilinks array of dimension M
%For M=2 c(1)=<k10>=<k01> c(2)=<k11>
%For M=2 c(1)=<k100>=<k010>=<k001> c(2)=<k110>=<k011>=<k101> c(3)=<k111>
%
% The output is  a cell array A of dimension M
% A{1} is the adjacency matrix of the first layer 
%          
% A{2} is the adjacency matrix of the second layer 
%
% A{3} is the adjacency matric of the third layer (for M=3) 
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
% [1] G. Bianconi
% Statistical mechanics of multiplex networks: Entropy and overlap." 
% Physical Review E 87, no. 6 (2013): 062806.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = Poisson_Multiplex_Model_Multilinks(N,M,c)

if (M==2)
    x=rand(N,N);
    y10=(ones(N,N)).*(x<(c(1)/N));
    y11=(x>(c(1)/N)).*(x<(c(1)+c(2)/N));
    y01=(x<(c(1)+c(2)/N)).*(x<(2*c(1)+c(2))/N):
    
    A{1}=y10+y11;
    A{2}=y01+y11;
end

if (M==3)
    
    x=rand(N,N);
    y100=(ones(N,N)).*(x<(c(1)/N));
    y110=(x>(c(1)/N)).*(x<(c(1)+c(2)/N));
    y010=(x<(c(1)+c(2)/N)).*(x<(2*c(1)+c(2))/N):
    y001=(x>(2*c(1)+c(2))/N).*(x<(3*c(1)+c(2))/N);
    y101=(x>(3*c(1)+c(2))/N).*(x<(3*c(1)+2*c(2))/N);
    y011=(x>(3*c(1)+2*c(2))/N).*(x<(3*c(1)+3*c(2))/N);
    y111=(x<(3*c(1)+3*c(2))/N).*(x<(3*c(1)+3*c(2)+c(3))/N);
    
    A{1}=y100+y110+y101+y111;
    A{2}=y010+y110+y011+y111;
    A{3}=y001+y101+y011+y111;
end
end



 
