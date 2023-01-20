
# check the quality of the data

docker run -it --rm -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized/All:/data:ro \
                    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/QC/Report_all:/out \
                    poldracklab/mriqc:latest /data /out participant

docker run -it --rm -v /Users/wang/Desktop/Research_projects/BBSC/Functional/Data/Reorganized/All:/data:ro \
                    -v /Users/wang/Desktop/Research_projects/BBSC/Functional/QC/Report_all:/out \
                   poldracklab/mriqc:latest /data /out group
