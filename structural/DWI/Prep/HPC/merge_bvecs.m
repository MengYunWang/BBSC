function merge_bvecs(fbval1,fbval2,fbvec1,fbvec2)
% Ivan I. Maximov, Ullevål
% 23.11.2018, ver 1.0

bl1 = load(fbval1);
bl2 = load(fbval2);
bv1 = load(fbvec1);
bv2 = load(fbvec2);

[fpath,fname,fext] = fileparts(fbval1);

bvec = [bv1 bv2];
bval = [bl1 bl2];

NN = size(bvec,2);

h = fopen([fpath '/data.bvec'],'w');
for i = 1:3
    for j = 1:NN
        fprintf(h,'%1.10f ',bvec(i,j));
    end
    fprintf(h,'\n');
end
fclose(h);

dlmwrite([fpath '/data.bval'],bval,'delimiter','\t','precision','%i');

end
