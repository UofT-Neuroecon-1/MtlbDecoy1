clear
load('simulout3.mat')
for s=1:40
    for t=1:250
        data(s).X{t}=V(t,:,s)';
    end
    data(s).y=choice{s};
    data(s).J=3;
    data(s).K=1;
end
save ExampleData3 'data' 'par'



clear
load('simuloutSS.mat')

for s=1:30
    for t=1:270
        X=V((s-1)*270+t,:)';
        X(isnan(X)) = [];
        data(s).X{t}=X;
        data(s).J(t)=length(X);
        data(s).K=1;
    end
    data(s).y=choice{1}((s-1)*270+1:s*270);  
end

load trueparSS.mat
par.scale=.5;
par.omega=omega; %Normalized?
par.sigma=sigma; %Gain control
par.a=a;
par.b=b; %power inside sum of denominator
par.kappa=1;

save ExampleDataSS 'data' 'par'