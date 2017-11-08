function gf = estimateGradient(f,x)
	epsilon = 10^(-3); % small epsilon as step for gradient estimation
	gf = zeros(numel(f(x)),numel(x));

	for i = 1:numel(x)
	    x_temp = x;
	    x_temp(i) = x_temp(i)+epsilon;
	    gf(:,i) = (f(x_temp)-f(x))/epsilon;
	end

	if numel(f(x))==1 % if f is real-valued function
	    gf = gf'; % should return a column vector
	end
end