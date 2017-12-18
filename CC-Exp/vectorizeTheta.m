function [ vecTheta ] = vectorizeTheta( Theta )
%VECTORIZETHETA Convert the cells of Theta into a G x P Matrix of vectors
%(a 3 dimensional matrix)
G = size(Theta,1);
P = size(Theta,2);
vecTheta = zeros(G,P,6);
for g=1:G
    for p=1:P
        theta = Theta{g,p};
        vecTheta(g,p,1:size(theta,2)) = theta;
    end
end

end

