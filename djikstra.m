%--------------------------------------------------------------------------
%* Function: Djikstra's Algorithm w/ Node Generation
%* Jeremy Kindseth
%* 12/21/2014
%* 
%* Notes: Input consists of total number of nodes (nPoints), the dimensions
%* of the space (xDim, yDim), a decay factor for determining nodal
%* connectivity (lambda ~0.05-0.5, lower is less connectivity), start node
%* (sNode) and finish node (fNode).  Output is optimal nodal path (path)
%*  and size of priority queue when algorithm terminates (pqSize)
%*
%--------------------------------------------------------------------------

function [path pqSize] = djikstra2(nPoints,xDim,yDim,lambda,sNode,fNode)

clc; clearvars -except nPoints xDim yDim lambda sNode fNode;
pts=nPoints;            
dim1=xDim;
dim2=yDim;
start=sNode;
last=fNode;

dim = 2;                                                            % Dimensionality
A = zeros(pts,dim); 

%% Build Node Network %%
i=0;
while(i<dim)
  i=i+1;
  s=strcat('dim',num2str(i));
  A(:,i)=eval(s).*rand(pts,1);
end  

G=zeros(pts);
dist=zeros(pts,pts);                                                % Distance matrix
output=0:1;

for i=1:pts                                                         % i = row
    j=0;
    while(j<pts)                                                    % j = column
        j=j+1;
        dist(i,j)=pdist([A(i,:);A(j,:)],'euclidean');
    end
end

sortDist=sort(dist(i,:));
maxDist=sortDist(pts);
minDist=sortDist(2);
diff=0.5*(maxDist-minDist);
lam=-1/(lambda*maxDist);                                            % exponential decay term
z=exp(lam.*dist(i,:));
    
k=0;
while(k<pts)
    k=k+1;
    for g=k:pts
        if (z(k)~=1)
           G(k,g)=randsample(output,1,true,[(1-z(k)) z(k)]);
        else
           G(i,k)=0;
        end
    end
    while(sum(G(k,:))==0)
        num=randi([k pts],1);                                       % If node is unconnected, give it one connection
        if(num~=k)
            G(k,num)=1;                                              
        elseif(k==pts)
            break            
        end
    end
end

for z=1:pts
    for o=1:pts
        G(o,z)=G(z,o);
    end
end

 %% Plot Nodes %%
 figure;
 axes;
 hold on;
 r=1.5;                                                         % Node size (radius)
 lbl = cellstr(num2str((1:pts)'));                              % Node labels
 title('Node Plot')
 gplot(G,A); hold on;
 xlabel('X Axis (Distance)');
 ylabel('Y Axis (Distance)');
 theta = linspace(0, 2*pi, 30)';
 x=A(:,1);
 y=A(:,2);
 xc = bsxfun(@plus, r .* cos(theta), x');
 yc = bsxfun(@plus, r .* sin(theta), y');
 patch(xc, yc, 'w');
 text(x,y,lbl);
 axis equal;

 %% Djikstra's Algorithm %%
 tic                                                            % Start clock to time execution
 cn=start;                                                      % Current node (set to start)
 cnindex=1;                                                     % Index of array of structs
 pq=struct('node',start,'visited',1,'value',0,'parent',start);  % Start Priority Queue
 % Find connected nodes
 while(1)
     for i=1:pts
         sz=length(pq);
         if(G(cn,i)==1) % Is node connected?
             flag=0;
             for j=1:sz % Does node exist?
                 if(pq(j).node==i && pq(j).visited==0)
                     flag=1;
                     if(pq(j).value>(pq(cnindex).value+dist(cn,i))) % If distance from current node is smaller than stored node value, update
                         pq(j).value=pq(cnindex).value+dist(cn,i);  % Update value
                         pq(j).parent=cn;
                     end
                 elseif(pq(j).node==i && pq(j).visited==1)
                     flag=1;
                 end
             end
             
             if(flag==0) % Node did not exist, create new node
                 pq(sz+1)=struct('node',i,'visited',0,'value',(pq(cnindex).value+dist(cn,i)),'parent',cn); % Create new structure
             end 
         end
     end
     % Find smallest value from node that hasn't been visited
     sz=length(pq);
     smallest=Inf;
     pqindex=Inf;
     for m=1:sz
         if(pq(m).visited==0 && pq(m).value<smallest) % If node has not been visited and value is smaller than stored smallest
             smallest=pq(m).value;
             cn=pq(m).node;  % Update current node
             cnindex=m;
         end
     end       
     % Mark new current node as visited
     pq(cnindex).visited=1;
     % Check to see if you've visited the goal node
     if(cn==last)
         % Find node path
         path=[cn];   % Initialize path array
         cnt=1;
         while(cn~=start)
             cnt=cnt+1;
             cn=pq(cnindex).parent;
             path=[cn path];  
             for y=1:length(pq)
                 if(pq(y).node==cn)
                     cnindex=y;
                     break
                 end
             end
         end
         toc                                                    % End timer
         pqSize=length(pq);                                     % Size of pq
         break                                                  % Full path found
     end
 end
end
 
 
     
 
 
 