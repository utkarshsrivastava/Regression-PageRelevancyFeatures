function [test_ans]= test_cfs(c,testlen,M,lambda,test_data,test_targetdata,mu,s_sq_inv)

phi_test=ones(testlen,c);
c=length(test_data);
phit_validate=ones(c,testlen);

	for i = 1:M-1
    for j =1:length(test_data)
	x(i)=test_data(j,:);
            phi(i,j)=exp(-0.5*transpose(x(i)-mu(j))*(s_sq_inv)*(x(i)-mu(j)));
        end
    end


phit_test=transpose(phi_test);
w=pinv(phit_test*phi_test+lambda*eye(c))*phit_test*test_targetdata;

Ermstest= sqrt(transpose(phi_test*w-test_targetdata)*(phi_test*w-test_targetdata)/testlen);
test_ans=Ermstest;

end
