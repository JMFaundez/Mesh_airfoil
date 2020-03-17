function [error, xu] = abs_error(x1, v1, x2, v2)
	v_int = interp1(x2,v2,x1);
	error = abs((v_int-v1)./v1)*100;
	xu = x1;
end
