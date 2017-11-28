function [ S ] = fsym( X )
%S = [X,X(:,end-1:-1:2); X(end-1:-1:2,:), X(end-1:-1:2,end-1:-1:2)];
S = [X,fliplr(X);flipud(X),fliplr(flipud(X))];
end

