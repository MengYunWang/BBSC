
# check the quality of the data

docker run -it --rm -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized/All:/data:ro \
                    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/QC/Report_all:/out \
                    nipreps/mriqc:23.1.0 /data /out participant

docker run -it --rm -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized/All:/data:ro \
                    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/QC/Report_all:/out \
                   nipreps/mriqc:23.1.0 /data /out group
