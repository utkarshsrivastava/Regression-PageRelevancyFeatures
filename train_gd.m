function[rms_gdtrain,rms_gdvalid,rms_gdtest]=train_gd(M)

load project1_data.mat
format long  
[rn,nc]=size(A);



trainlen = floor(0.8*rn); 
validlen=floor(0.1*rn); 
testlen=floor(0.1*rn);

%Seperate out the target vectors
validation_targetdata=A(trainlen+1:trainlen+1+validlen,1);
test_targetdata=A(trainlen+validlen+1:trainlen+validlen+testlen,1);
train_targetdata=A(1:trainlen,1); % Store the target vector for training data

train_data = A(1:trainlen,2:47);  
validate_data= A(trainlen+1:trainlen+1+validlen,2:47); 
test_data=A(trainlen+validlen+1:trainlen+validlen+testlen,2:47); % Load the test data

eta=0.5;

mumat=zeroes(F,M);
x=zeroes(F,1);
rndmu=M;
for i = 1: M
	iter=randi(N)
	for j=1:F
	mu(i)=IPData(iter);
	end
end
 


%calculate phi
rndmu=M;

s_sq_inv=inverse(s*s*(eye(N))));

phi(:,1)=1;
 for(i=1:N)
  
  for (j=2:M)
	
  x(i,:)=transpose(IPData(1,i));
  phi(i,j)=exp(-0.5*transpose(x(i)-mu(j))*(s_sq_inv)*(x(i)-mu(j)));
	end
  end

save mu_cfs mu;
save s_gd s;

%Training 
M=10;
features=(M-1)*46+1;
num_samples=length(train_data);
w=zeros(M,1);

for i=1:M-1         
w(i)=rndi(1);% phi matrix how ?
end


for i=1:M-1         
 w(i+1)=w(i)+eta*(train_targetdata(i)-(transpose(w(i))*(phi(i,:))))*transpose(phi(i,:));
end
count=0;
Ermsmat=zeros(100,1);
Ermstrain= sqrt(transpose(phi*w-train_targetdata)*(phi*w-train_targetdata)/trainlen);
Ermstrain1=0;
errorsg=abs(Ermstrain1-Ermstrain);
errorsg1=errorsg;
while count<100 or (errorsg1< 0.1 do
Ermstrain=Ermstrain1;
Ermstrain1= sqrt(transpose(phi*w-train_targetdata)*(phi*w-train_targetdata)/trainlen);
count=count+1;
Ermsmat(count,1)=Ermstrain1;
errorsg=errorsg1;
errorsg1=abs(Ermstrain1-Ermstrain);

if (errorsg1>errorsg) 
then
eta=eta*0.5;
end_if

end_while
%Ermstrain= sqrt(transpose(phi*w-train_targetdata)*(phi*w-train_targetdata)/trainlen);
rms_gdtrain(l)=Ermstrain;
save w_cfs w;
  
%Validation
for i=1:features-1         
        w(i+1)=w(i)+ eta*(validation_targetdata(i)-(transpose(w(i))*(phi(i,:))))*transpose(phi(i,:));
end

phit=phi(1:validlen,:);
Ermsvalid= sqrt(transpose(phit*w- validation_targetdata(1:validlen,:))*(phit*w-validation_targetdata(1:validlen,:))/validlen);
rms_gdvalidation(l)=Ermsvalid;

rms_gdtest= test_gd(testlen,M,test_data,test_targetdata,coeffs,l,features,eta);
end


%end


