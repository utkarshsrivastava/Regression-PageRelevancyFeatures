function [testans_gd]= test_gd(testlen,M,test_data,test_targetdata,mu,s_sq_inv,l,F,eta)
test_phi=ones(testlen,46);
trans_test_phi=ones(46,testlen);
w=zeros(F,1);
for i = 1:M-1
    for j =1:length(test_data)
	x(i)=test_data(j,:);
            test_phi(i,j) =exp(-0.5*transpose(x(i)-mu(j))*(s_sq_inv)*(x(i)-mu(j));
        end
    end
end

for i=1:F-1         
        w(i+1)=w(i)+ eta*(test_targetdata(i)-(transpose(w(i))*(test_phi(i,:))))*transpose(test_phi(i,:));
end
phit=test_phi(1:testlen);
Ermstest= sqrt(transpose(w*phit- test_targetdata(1:testlen))*(phit*w-test_targetdata(1:testlen))/testlen);
testans_gd=Ermstest;
end