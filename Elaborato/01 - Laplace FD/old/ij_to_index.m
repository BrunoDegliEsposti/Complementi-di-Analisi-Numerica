function [index] = ij_to_index(i,j,Ntheta)
    if i == 0
        i = Ntheta;
    elseif i == Ntheta+1
        i = 1;
    end
    index = i+(j-2)*Ntheta;
end

