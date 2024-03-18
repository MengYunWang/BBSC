# make command list
# then make nii list
if (!require("pacman")) install.packages("pacman")
pacman::p_load(sqldf)
# 1) prep commands for dcm to nii conversion (requiring only a varying input file)
allfiles = list.dirs(path = "/cluster/projects/p33/users/maxk/BBSC/structural", full.names = T)
filesvector = data.frame(allfiles)
dwi_files = sqldf("select * from filesvector where allfiles LIKE '%DTI%'")
cmd_list = c()
for (i in 1:nrow(dwi_files)){
  cmd_list[i] = paste("sh build_nii2.sh ", dwi_files$allfiles[i], sep = "")
}
write.table(cmd_list,"/cluster/projects/p33/users/maxk/BBSC/cmd_list1.txt", row.names = F, quote = F, col.names = F)
#
# this command list was used to create one .nii per session (per participant)
# the locations of the resulting .nii files are listed in the following list
#
# 2) now, make a list of the resulting nifti files
file_names = list()
for (list in 1:nrow(dwi_files)){
  file_names[[list]] = list.files(path = paste(dwi_files$allfiles[list]), pattern = "nii", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/nii_list1.txt", row.names = F, quote = F, col.names = F)
# and for each subject separately
## sub1
sub1files = list.dirs(path = "/cluster/projects/p33/users/maxk/BBSC/structural/sub1/", full.names = T)
sub2files = list.dirs(path = "/cluster/projects/p33/users/maxk/BBSC/structural/sub2/", full.names = T)
sub3files = list.dirs(path = "/cluster/projects/p33/users/maxk/BBSC/structural/sub3/", full.names = T)
file_names = list()
for (list in 1:length(sub1files)){
  file_names[[list]] = list.files(path = paste(sub1files[list]), pattern = "nii", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub1_nii_list1.txt", row.names = F, quote = F, col.names = F)
## sub2
file_names = list()
for (list in 1:length(sub2files)){
  file_names[[list]] = list.files(path = paste(sub2files[list]), pattern = "nii", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub2_nii_list1.txt", row.names = F, quote = F, col.names = F)
## sub3
file_names = list()
for (list in 1:length(sub3files)){
  file_names[[list]] = list.files(path = paste(sub3files[list]), pattern = "nii", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub3_nii_list1.txt", row.names = F, quote = F, col.names = F)

# make a list of the b-val files
file_names = list()
for (list in 1:length(sub1files)){
  file_names[[list]] = list.files(path = paste(sub1files[list]), pattern = "bval", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub1_bval_list1.txt", row.names = F, quote = F, col.names = F)
## sub2
file_names = list()
for (list in 1:length(sub2files)){
  file_names[[list]] = list.files(path = paste(sub2files[list]), pattern = "bval", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub2_bval_list1.txt", row.names = F, quote = F, col.names = F)
## sub3
file_names = list()
for (list in 1:length(sub3files)){
  file_names[[list]] = list.files(path = paste(sub3files[list]), pattern = "bval", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub3_bval_list1.txt", row.names = F, quote = F, col.names = F)



# make a list of the bvec files
file_names = list()
for (list in 1:length(sub1files)){
  file_names[[list]] = list.files(path = paste(sub1files[list]), pattern = "bvec", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub1_bvec_list1.txt", row.names = F, quote = F, col.names = F)
## sub2
file_names = list()
for (list in 1:length(sub2files)){
  file_names[[list]] = list.files(path = paste(sub2files[list]), pattern = "bvec", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub2_bvec_list1.txt", row.names = F, quote = F, col.names = F)
## sub3
file_names = list()
for (list in 1:length(sub3files)){
  file_names[[list]] = list.files(path = paste(sub3files[list]), pattern = "bvec", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub3_bvec_list1.txt", row.names = F, quote = F, col.names = F)


# make a list of the json files
file_names = list()
for (list in 1:length(sub1files)){
  file_names[[list]] = list.files(path = paste(sub1files[list]), pattern = "json", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub1_json_list1.txt", row.names = F, quote = F, col.names = F)
## sub2
file_names = list()
for (list in 1:length(sub2files)){
  file_names[[list]] = list.files(path = paste(sub2files[list]), pattern = "json", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub2_json_list1.txt", row.names = F, quote = F, col.names = F)
## sub3
file_names = list()
for (list in 1:length(sub3files)){
  file_names[[list]] = list.files(path = paste(sub3files[list]), pattern = "json", full.names = T)
}
file_names = data.frame(do.call(rbind,file_names))
write.table(file_names,"/cluster/projects/p33/users/maxk/BBSC/sub3_json_list1.txt", row.names = F, quote = F, col.names = F)


#
#
#
#
#
# potentially duplicated nifti, a.bval, a.bvec, a.json files which get a letter suffix can be deleted
#
# make list of .nii files which can be deleted
deleteFiles = data.frame(do.call(rbind,file_names))
delete_files = sqldf("select * from deleteFiles where X2 LIKE '%a.nii%'")
del = c()
for (i in 1:nrow(delete_files)){
  del[i] = (paste("rm",delete_files$X2[i]))
}
write.table(del,"/cluster/projects/p33/users/maxk/BBSC/delete_list1.txt", row.names = F, quote = F, col.names = F)

# bval files that can be deleted
file_names = c()
for (list in 1:nrow(dwi_files)){
  file_names[list] = list.files(path = paste(dwi_files$allfiles[list]), pattern = "a.bval", full.names = T)
}
file_names = data.frame(file_names)
delete_files = sqldf("select * from file_names where file_names LIKE '%a.bval%'")
del = c()
for (i in 1:nrow(delete_files)){
  del[i] = (paste("rm",delete_files$file_names[i]))
}
write.table(del,"/cluster/projects/p33/users/maxk/BBSC/delete_list2.txt", row.names = F, quote = F, col.names = F)

# bvec files that can be deleted
file_names = c()
for (list in 1:nrow(dwi_files)){
  file_names[list] = list.files(path = paste(dwi_files$allfiles[list]), pattern = "a.bvec", full.names = T)
}
file_names = data.frame(file_names)
delete_files = sqldf("select * from file_names where file_names LIKE '%a.bvec%'")
del = c()
for (i in 1:nrow(delete_files)){
  del[i] = (paste("rm",delete_files$file_names[i]))
}
write.table(del,"/cluster/projects/p33/users/maxk/BBSC/delete_list3.txt", row.names = F, quote = F, col.names = F)


# json files that can be deleted
file_names = c()
for (list in 1:nrow(dwi_files)){
  file_names[list] = list.files(path = paste(dwi_files$allfiles[list]), pattern = "a.json", full.names = T)
}
file_names = data.frame(file_names)
delete_files = sqldf("select * from file_names where file_names LIKE '%a.json%'")
del = c()
for (i in 1:nrow(delete_files)){
  del[i] = (paste("rm",delete_files$file_names[i]))
}
write.table(del,"/cluster/projects/p33/users/maxk/BBSC/delete_list4.txt", row.names = F, quote = F, col.names = F)


print("Script finished")

