x=0.1:0.025:0.9;
allcvibint=0;
for j=1:900
    for i=1:30
        A=allcvibampl(i,j);
        FWHM=allcvibfwhm(i,j);
        Peak=allcvibpeaks(i,j);
    allcvibint(i,j)=sum(gaussian(650:1:820,allcvibampl(i,j),allcvibfwhm(i,j),allcvibpeaks(i,j)));
    %plot(gaussian(650:1:820,allcvibampl(i,j),allcvibfwhm(i,j),allcvibpeaks(i,j)))
    Int=allcvibint(i,j);
    end
    
end
fractions=allcvibint(goodvibs)./(allcint1(goodvibs)+allcint2(goodvibs));
histogrammutant_2=histc(fractions,x);
figure(30)
bar(x, histogrammutant_2)