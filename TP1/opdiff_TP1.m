%%

gradx=@(x) [zeros(size(x,1),1),x(:,2:end)-x(:,1:end-1)];
grady=@(x) [zeros(1,size(x,2));x(2:end,:)-x(1:end-1,:)];

prepdivx=@(x) [x(:,2) , x(:,3:end)-x(:,2:end-1),-x(:,end)];
prepdivy=@(x) [x(2,:) ; x(3:end,:)-x(2:end-1,:);-x(end,:)];

div=@(x,y) prepdivx(x)+prepdivy(y);

