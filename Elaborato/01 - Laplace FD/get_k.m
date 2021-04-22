function [k] = get_k(i,j,Ntheta)
    % il blocco if serve a rendere ciclici gli indici j. questo tornerà
    % più utile avanti per imporre le condizioni al bordo periodiche.
    if (j == Ntheta+1)
        j = 1;
    elseif (j == 0)
        j = Ntheta;
    end
    k = (i-1)*Ntheta + j;
end
