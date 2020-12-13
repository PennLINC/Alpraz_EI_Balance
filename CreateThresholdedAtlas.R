library(RNifti)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  print("missing argumens (atlas, threshold)")
} else {
  mask_file <- args[1]
  threshold <- as.numeric(args[2])
  }

#Read in the ROI data
atlas_means<-read.table('/cbica/projects/alpraz_EI/input/atlases/atlas_means.txt',header = F,stringsAsFactors = F,quote = "")
atlas_means<-atlas_means[,2:dim(atlas_means)[2]]
atlas_means<-t(atlas_means)
atlas_means[,1] <- sub(x=atlas_means[,1],pattern = "Mean_",replacement = "")
atlas_means<-as.data.frame(atlas_means,row.names = NULL)
colnames(atlas_means)=c("index","value")
rownames(atlas_means)=NULL
atlas_means$index<-as.numeric(levels(atlas_means$index))[atlas_means$index]
atlas_means$value<-as.numeric(levels(atlas_means$value))[atlas_means$value]

# Find atlas regions that do not meet the threshold
bad_rois <- atlas_means$index[atlas_means$value<threshold]

# Load the atlas
atlas <- readNifti(mask_file)
# Make a new image to store output
output <- atlas

# Zero out all atlas regions that do not meet threshold
output[output %in% bad_rois] = 0

# write output
writeNifti(output,sprintf('/cbica/projects/alpraz_EI/input/atlases/atlas_threshold_%1.2f.nii.gz',threshold))
