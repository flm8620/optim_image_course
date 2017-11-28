function [ S ] = fsym2( X )
S = [X,X(:,end-1:-1:2); X(end-1:-1:2,:), X(end-1:-1:2,end-1:-1:2)];
end

