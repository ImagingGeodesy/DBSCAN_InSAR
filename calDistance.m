%% calculate distance between points in the matrix
% m points, n dimensions
function dis = calDistance(x)
    [m,~] = size(x);
    dis = zeros(m,m);
    
    for i = 1:m
        for j = i:m
            dis(i,j) = sqrt(sum((x(i,:)-x(j,:)).^2));
            dis(j,i) = dis(i,j);
        end
    end
end