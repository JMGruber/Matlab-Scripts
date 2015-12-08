readdir='O:\Michael\2015\TOM\LHCII homogeneity experiment\Sample TOM';
writedir='O:\Michael\2015\TOM\LHCII homogeneity experiment\LHCII CP24 KO control\spectral analysis';
%bg=dlmread(fullfile(readdir,'bgall.txt'));

includefiles=[61 127 211 261];
% WT includefiles=[79 80 126 142 160 164 169 168 200 231]

%Meanval=zeros(225,1);
Meanval=dlmread(fullfile('O:\Michael\2015\TOM\LHCII homogeneity experiment\LHCII CP24 KO control',['bgfile.txt']));
% for specnumber=1:300;
% j=1;
%     while (j<=length(includefiles))&&(includefiles(j)<=specnumber)
%         if specnumber==includefiles(j)
%             mat=dlmread(fullfile('O:\Michael\2015\TOM\LHCII homogeneity experiment\Sample TOM',['spec' int2str(specnumber)]));
%             len=length(mat);
%             Meanval=zeros(len,1);
%             Meanval=(Meanval+mean(mat(:,2:size(mat,2)),2))/2;
%         end
%         j=j+1;
%     end
%     
% 
% 
% end
j
mat_avg1=zeros(len,1);
mat_avg2=zeros(len,1);


for specnumber=37;
    mat=dlmread(fullfile(readdir,['spec' int2str(specnumber)]));
    
    mat_2=zeros(len,1);
    mat_3=zeros(len,1);
    Length2=length(mat(1,:));
    for j=2:Length2;
        Sum1=sum(mat(:,j)-Meanval);
        Sum2=sum(mat(1:45,j)-Meanval(1:45));
        if Sum1>300;
            mat_2=mat_2+mat(:,j)-Meanval;
        if Sum2<50;
        %if sum(mat(:,j)-Meanval)>1000;
            mat_3=mat_3+mat(:,j)-Meanval;
        end
        %mat_2=smooth(mat_2,3);
        end
    end
    %plot(mat(:,1),mat_2)
    mat_avg1=mat_avg1+mat_2;
    mat_avg2=mat_avg2+mat_3;
    
    
   
 
   % h=image(mat(:,1),1:size(mat,2),mat(:,2:end)');
   % axis([600 800 0.5 size(mat,2)+.5])
   % xlabel('Wavelength (nm)');
   % ylabel('Illumination time (s)');
   % colormap(jet(min(round(max(max(mat(:,2:end)))*5/6),256)));
   % saveas(h,fullfile(writedir,['contourplot' int2str(specnumber) '.jpg']));
end
%mat_avg=smooth(mat_avg,15);
%mat_avg1=mat_avg1/norm(mat_avg1);
%mat_avg2=mat_avg2/norm(mat_avg2);
figure(1)
plot(mat(:,1),mat_avg1,mat(:,1),mat_avg2)
                            
   