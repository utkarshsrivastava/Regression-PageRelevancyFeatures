
function[rms_train,rms_valid,rms_test]=  train_cfs(M,lambda)
load project1_data.mat
%Linear Regression with basis function
format long  %Format all the computation in long
[rn,nc]=size(A);
%Init
%Number of features is 46

%Split Dataset into training, validation and testing set
trainlen = floor(0.8*rn); % 80% Training data
validlen=floor(0.1*rn); %10% Validation data
testlen=floor(0.1*rn);
%Seperate out the target vectors
train_targetdata=A(1:trainlen,1); % Store the target vector for training data
validation_targetdata=A(trainlen+1:trainlen+1+validlen,1);
test_targetdata=A(trainlen+validlen+1:trainlen+validlen+testlen,1);


train_data = A(1:trainlen,2:47);  
validate_data= A(trainlen+1:trainlen+1+validlen,2:47); 
test_data=A(trainlen+validlen+1:trainlen+validlen+testlen,2:47);

s=0.5;
s_sq_inv=inverse(s*s*(eye(N))));
mumat=zeroes(F,M);
x=zeroes(F,1);
rndmu=M;
%mumat=Util.getRandomMu();
for i = 1: M
	iter=randi(N);
	for j=1:F
	mu(i,j)=IPData(iter);
	end
	end
mumat=transpose(mu);
save('mu_cfs',mumat);
save('s_cfs',s);

phi=zeroes(M,N);
phi(:,1)=1;
 for(i=1:N)
  k=trainlen/i;
  for (j=2:M)
	
  x(i,:)=transpose(IPData(1,i));
  phi(i,j)=exp(-0.5*transpose(x(i)-mu(j))*(s_sq_inv)*(x(i)-mu(j)));
	end
  end

%w=calculate w(M,N,phi,lambda,target)  
  w=zeroes(M,1);
  tphi=zeroes(M,N);
  tphi=transpose(phi);
  lambI=lambda*eye(M);
  w=pinv(tphi*phi+lambI)*tphi*target;
Erms=0.0;
Erms= sqrt(transpose(phi*w-train_targetdata)*(phi*w-train_targetdata)/trainlen);
%rms_train=Ermstrain;

	%Validation
val_phi=ones(validlen,M);
transval_phi=ones(M,validlen);
for i = 1:M-1
    for j =1:length(validate_data)
            phi(i,j)=exp(-0.5*transpose(x(i)-mu(j))*(s_sq_inv)*(x(i)-mu(j)));
        end
    end


transval_phi=transpose(val_phi);
w=pinv(transval_phi*val_phi+lambda*eye(c))*transval_phi*validation_targetdata;

Ermsval= sqrt(transpose(val_phi*w-validation_targetdata)*(val_phi*w-validation_targetdata)/validlen);
rms_valid=Ermsval;
rms_test= test_cfs(c,testlen,M,lambda,test_data,test_targetdata,mumat,s);
end