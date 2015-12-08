function[B]=imageshift(C);
B=C(1:101,1:100);
%B = imrotate(B,90);
for i=1:2:100;
    B(i,:)=circshift(B(i,:),[0 2]);
end;
B = flipud(B);
%B = imrotate(B,90);
figure(1)
image(B, 'CDataMapping', 'scaled')
