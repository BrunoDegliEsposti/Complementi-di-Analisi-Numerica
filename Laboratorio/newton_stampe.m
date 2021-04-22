function [x,fx,k,flag] = newton_stampe(x,f,J,tolla,tollr,tollf,kmax,i_stampe)
%NEWTON_STAMPE Metodo di Newton per sistemi non lineari
    n = length(x);
	fx = f(x);
	if norm(fx) <= tollf
		flag = 1;
        k = 0;
		return
	end
	dxold = ones(n, 1);
    if i_stampe == 1
        fprintf(['\n k   xk(1)\t\t    norm(f(xk))\t\t    diff', ...
                 '\t     rapp\t rapp^2\n']);
    end
	for k = 1:kmax
		Jx = J(x);
		if any(any(isinf(J(x)))) || rcond(Jx) < eps
			flag = -1;  % La matrice jacobiana è mal condizionata
			return
		end
		dx = Jx \ fx;
        x = x - dx;
        fx = f(x);
        if i_stampe == 1
			fprintf(' %02u %22.15e %22.15e %16.8e %11.3e %11.3e\n', ...
                    k, x(1), norm(fx), norm(dx), norm(dx)/norm(dxold), ...
                    norm(dx)/(norm(dxold)^2));
        end
        if norm(fx) <= tollf
            flag = 2;
            return
        end
		if norm(dx) <= tolla + tollr*norm(x)
			flag = 1;
			return
		end
		dxold = dx;
	end
	flag = 0;
end

