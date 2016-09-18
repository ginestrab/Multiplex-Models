% %%%%%%%%%%%%%Scale-free Multiplex Model Multilinks %%%%%%%%%%%%%%%%%%%
% This code generates a  undirected unweighted duplex or triplex network with 
% Poisson distribution of multilinks
%
% INPUTS: 
%
% N total number of nodes in the multiplex
% M total number of layers (M=2 or M=3)
% gamma power-law exponent of multilinks: array of dimension M
%For M=2 
%    -gamma(1) power-law exponent multidegrees k10, k01
%    -gamma(2) power-law exponent multidegrees k11
%For M=3 
%    -gamma(1) power-law exponent multidgrees k100, k010, k001
%    -gamma(1) power-law exponent multidgrees k110, k011, k101
%    -gamma(3) power-law exponent multidgree k111
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

function [A] = Scalefree_Multiplex_Model_Multilinks(N,M,gamma)

if (M==2),
    for i=1:N,
        k10(i)=rand(1).^(-1/(gamma(1)-1));
        
        while(k10(i)>sqrt(N))
            k10(i)=rand(1).^(-1/(gamma(1)-1));
        end
        k01(i)=rand(1).^(-1/(gamma(1)-1));
        while(k01(i)>sqrt(N))
            k01(i)=rand(1).^(-1/(gamma(1)-1));
        end
        k11(i)=rand(1).^(-1/(gamma(2)-1));
        while(k11(i)>sqrt(N))
            k11(i)=rand(1).^(-1/(gamma(2)-1));
        end
    end
    k01=k01';
    k10=k10';
    k11=k11';

    x=tril(rand(N,N));
    y10=tril((ones(N,N)).*(x<(k10*k10'/sum(k10))));
    y11=tril((x>(k10*k10'/sum(k10))).*(x<((k10*k10'/sum(k10))+(k11*k11'/sum(k11)))));
    y01=tril((x>((k10*k10'/sum(k10))+(k11*k11'/sum(k11)))).*(x<((k10*k10'/sum(k10))+(k11*k11'/sum(k11))+(k01*k01'/sum(k01)))));
        
    A{1}=y10+y11;
    A{2}=y01+y11;
    A{1}=A{1}+A{1}';
    A{2}=A{2}+A{2}';
end

if (M==3),

        for i=1:N,
        k100(i)=rand(1).^(-1/(gamma(1)-1));
        
        while(k100(i)>sqrt(N))
            k100(i)=rand(1).^(-1/(gamma(1)-1));
        end
        k010(i)=rand(1).^(-1/(gamma(1)-1));
        while(k010(i)>sqrt(N))
            k010(i)=rand(1).^(-1/(gamma(1)-1));
        end
        k110(i)=rand(1).^(-1/(gamma(2)-1));
        while(k110(i)>sqrt(N))
            k110(i)=rand(1).^(-1/(gamma(2)-1));
        end
          k101(i)=rand(1).^(-1/(gamma(2)-1));
        
        while(k101(i)>sqrt(N))
            k101(i)=rand(1).^(-1/(gamma(2)-1));
        end
        k011(i)=rand(1).^(-1/(gamma(2)-1));
        while(k011(i)>sqrt(N))
            k011(i)=rand(1).^(-1/(gamma(2)-1));
        end
        k111(i)=rand(1).^(-1/(gamma(3)-1));
        while(k111(i)>sqrt(N))
            k111(i)=rand(1).^(-1/(gamma(3)-1));
        end
        k001(i)=rand(1).^(-1/(gamma(1)-1));
        while(k001(i)>sqrt(N))
            k001(i)=rand(1).^(-1/(gamma(3)-1));
        end
    end
    k010=k010';
    k100=k100';
    k110=k110';
    k011=k011';
    k101=k101';
    k111=k111';
    k001=k001';

    p010=k010*k010'/sum(k010);
    p100=k100*k100'/sum(k100);
    p110=k110*k110'/sum(k110);
    p011=k011*k011'/sum(k011);
    p101=k101*k101'/sum(k101);
    p111=k111*k111'/sum(k111);
    p001=k001*k001'/sum(k001);
    
    x=tril(rand(N,N));
    y100=tril((ones(N,N)).*(x<p100));
    y110=tril((x>p100).*(x<(p100+p110)));
    y010=tril((x>(p100+p110)).*(x<(p100+p110+p010)));
    y001=tril((x>(p100+p110+p010)).*(x<(p100+p110+p010+p001)));
    y101=tril((x>(p100+p110+p010+p001)).*(x<(p100+p110+p010+p001+p101)));
    y011=tril((x>(p100+p110+p010+p001+p101)).*(x<(p100+p110+p010+p001+p101+p011)));
    y111=tril((x>(p100+p110+p010+p001+p101+p011)).*(x<(p100+p110+p010+p001+p101+p011+p111)));
    
    A{1}=y100+y110+y101+y111;
    A{2}=y010+y110+y011+y111;
    A{3}=y001+y101+y011+y111;
    A{1}=A{1}+A{1}';
    A{2}=A{2}+A{2}';
    A{3}=A{3}+A{3}';
end
end



 
