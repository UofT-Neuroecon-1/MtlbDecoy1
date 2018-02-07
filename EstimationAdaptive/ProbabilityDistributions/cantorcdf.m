
% v=rand(1,100);
% v=linspace(0,1,10001);
% 
% %Convert to base 3
% temp=v(:);
% N=zeros(size(temp));
% sz=size(v);
% ind=1;
% n=floor(log(1/eps)/log(3));
% clear M
% for i=1:n
%     temp2=temp.*3;
%     M(:,i)=floor(temp2);
%     temp=temp2-M(:,i);
% end
% 
% %Find first 1's in each row 
% [ignore,ind]=max(M'==1);
% 
% %Erase all values afterwards
% for i=1:length(temp)
%     M(i,ind(i)+1:end)=0;
% end
% 
% %Replace all 2 with 1
% M(M==2)=1;
% 
% %Convert from decimal
% M=M*(2.^(-(1:n)))';
% 
% plot(M,'.')


x=rand(10,1);  %numbers (0,1)
%Multiply to make integer
v=round(x*(2^52));
%Express in base 3
M=dec2base(v,3)
