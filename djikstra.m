clc; clear all;

% The data

pts = 800;              % Number of nodes
dim = 2;                % Dimensionality
A = zeros(pts,dim);     % Initialize
dim1 = 100;           % Linear space dimension
dim2 = 100;           % Linear space dimension

% Distribute points throughout space
i=0;
while(i<dim)
  i=i+1;
  s=strcat('dim',num2str(i));
  A(:,i)=eval(s).*rand(pts,1);
end  

G=zeros(pts);
dist=zeros(pts,pts);        % Distance matrix
output=0:1;

for i=1:pts                 % i = row
    j=0;
    while(j<pts)            % j = column
        j=j+1;
        dist(i,j)=pdist([A(i,:);A(j,:)],'euclidean');
    end
end

sortDist=sort(dist(i,:));
maxDist=sortDist(pts);
minDist=sortDist(2);
diff=0.5*(maxDist-minDist);
lam=-1/(0.1*maxDist);         % exponential decay term
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
        num=randi([k pts],1);
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



 x = [1 2 2 2 ]; % Node x coordinates
 y = [3 3 2 1];  % Node y coordinates
 nx = length(x);
 lbl = cellstr(num2str((1:pts)')); % Node labels
 r = 1.5; % radius of nodes
 adj = zeros(nx); % Adjacency matrix for edges,
 adj(1,2) = 1;    % values specify line width
 adj(1,3) = 1;
 adj(1,4) = 1;
 % The plot
 figure;
 axes;
 hold on;
 linewidth = unique(adj(adj>0));
 gplot(G,A); hold on;
 theta = linspace(0, 2*pi, 30)';
 x=A(:,1);
 y=A(:,2);
 xc = bsxfun(@plus, r .* cos(theta), x');
 %xc = bsxfun(@plus, r .* cos(theta), A(:,1));
 yc = bsxfun(@plus, r .* sin(theta), y');
 %yc = bsxfun(@plus, r .* sin(theta), A(:,2));
 patch(xc, yc, 'w');
 text(x,y,lbl);
 axis equal;
 
 
 tic
 % Djikstra's Algorithm
 % define start and end nodes
 start = 1;
 last = pts;
 cn=start;              % current node (set to start)
 cnindex=1;
 %pq=[linspace(1,pts,pts)' zeros(pts,1) inf(pts,1) start.*ones(pts,1)];   % Priority Queue: [(int)node index; (bool)visited; (int)node value; (string)path]
 %pq(start,2)=1;         % set start node as visited
 
 pq=struct('node',start,'visited',1,'value',0,'parent',start);  % Visited because is current node
 
 
 % Find connected nodes
 while(1)
     for i=1:pts
         sz=length(pq);
         if(G(cn,i)==1) % Is node connected?
             flag=0;
             for j=1:sz % Does node exist?
                 if(pq(j).node==i && pq(j).visited==0)
                     % Exists in PQ
                     flag=1;
                     if(pq(j).value>(pq(cnindex).value+dist(cn,i))) % If distance from current node is smaller than stored node value, update
                         pq(j).value=pq(cnindex).value+dist(cn,i); % Update value
                         pq(j).parent=cn;
                     end
                 elseif(pq(j).node==i && pq(j).visited==1)
                     flag=1;
                 end
             end
             
             if(flag==0) % Node did not exist, create new node
                 %pq(sz+1)=struct('node',i,'visited',[],'value',(pq(cnindex).value+dist(cn,i)),'parent',cn);   % Create new structure
                 pq(sz+1)=struct('node',i,'visited',0,'value',(pq(cnindex).value+dist(cn,i)),'parent',cn);
                 
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
         %pq(cnindex).value
         path
         toc
         break
     end

     %break
 end