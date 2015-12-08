readdir='J:\2014\Joshua and Herman\pH dependence PCD PCA\pH 5\';

name='refocus';
Refocus=[];
for i=1:1000
    if exist([readdir name int2str(i)],'file')
        Refocus=[Refocus i];
    end
end