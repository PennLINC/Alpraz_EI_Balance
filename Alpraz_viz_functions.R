# Visualization functions for alpraz prediction
# library(tidyverse)
# library(oro.nifti)
library(neurobase)
library(cowplot)
library(ggridges)
library(pracma)
# library(purrr)
# library(mgcv)
library(gratia)
# library(broom)
library(scales)
# library(kableExtra)
# library(RColorBrewer)
source('tools/RainCloudPlots/tutorial_R/R_rainclouds.R')
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
theme_replace(axis.text=element_text(colour = "black",size = font_size))
line_size <- 1.5
point_size <- 2

perm_communities <- function(permW,atlasname){
  # This function calculates the average W for each community.
  # First, we need to match features with their community assignments
  ## Load an example correlation matrix (names are node Idx)
  fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK_GSR/13783/1986/%s_CC_GSR_000.netcc",atlasname)
  templateCC <- as.matrix(read.table(fname,header = T,skip = 4))
  permW <- as.data.frame(abs(permW))
  perm_rankings <- mclapply(mc.cores = 8,X = permW,feat2mat,
                            atlasname=atlasname,
                            community_summary = T,
                            plots = F,
                            templateCC = templateCC,
                            returnRankingsOnly = T)
  
  community_summary <- perm_rankings
  return(community_summary)
}


feat2mat <- function(W,atlasname,community_summary = F,plots = T,templateCC = NULL,returnRankingsOnly = F,GABA=T,transmodality=T,surface=F){
  if (is.null(templateCC)) {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK_GSR/13783/1986/%s_CC_GSR_000.netcc",atlasname)
    templateCC <- as.matrix(read.table(fname,header = T,skip = 4))
  }
  
  #load template matrix (picked 13783 1986 at random)
  
  templateCC_vector = templateCC[upper.tri(templateCC,diag = F)]
  if (length(templateCC_vector == length(W))) {
    templateCC_vector = W
  } else {
    stop("Length of W does not match the supplied atlas")
  }
  newCC <- templateCC
  newCC[upper.tri(newCC,diag = T)] = 0
  newCC[upper.tri(newCC,diag = F)] = templateCC_vector
  newCC[lower.tri(newCC)] <- t(newCC)[lower.tri(newCC)]
  
 
  # Get labels
  
  labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_labels.csv',atlasname),header=F)
  if (atlasname %in% c("schaefer400x7_aal","schaefer200x7_aal","gordon333_aal","glasser360_aal")) {
    community_labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_node_labels.csv',
                                         gsub(x = atlasname,pattern = "_aal",replacement = ""))
                                 ,header=T)
    
    mat_labels <- data.frame(nums = colnames(templateCC))
    mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
    mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
    mat_labels <- mat_labels %>% left_join(community_labels,by = c("nums"="index"))
    mat_labels <- mat_labels %>%
      mutate(community_name = as.character(community_name),V2 = as.character(V2),node_names=as.character(node_names))
    mat_labels$community_name[is.na(mat_labels$community_name)]="AAL_subcortex"
    mat_labels$node_names[is.na(mat_labels$node_names)]= mat_labels$V2[is.na(mat_labels$node_names)]
    
    # Set names to node names
    colnames(newCC)=mat_labels$node_names
    rownames(newCC)=mat_labels$node_names
    #write to brain
    # mat2brain(newCC,atlasname,sprintf('%_weights.nii.gz',atlasname))
    
    # Arrange by community membership
    orderCC <- newCC[order(mat_labels$node_community),order(mat_labels$node_community)]
    newCC <- orderCC
    
    if (community_summary == T) {
      commCC <- newCC
      df_CC <- as.data.frame(commCC)
      rvals_abs <- rowSums(abs(commCC)) #get the sum of absolute values
      rvals <- rowSums(commCC) #get the sum of signed values
      rval_df <- data.frame(rvals = rvals,rvals_abs = rvals_abs, names = as.character(names(rvals)))
      rval_df <- rval_df %>% 
        left_join(mat_labels,by = c("names" = "node_names"))
      if (surface==T) {
        vec2surface(schaefer_df = rval_df,val_name = "rvals_abs",file_prefix = sprintf("feature_weights_abssum_%s",atlasname))
        #write files for brainsmash spatial correlation as well:
        for_brainsmash <- rval_df %>% arrange(nums) %>% select(rvals)
        write.table(for_brainsmash,"/cbica/projects/alpraz_EI/output/brainsmash/input_files/schaefer400x7_aal_threshold_0.95_SVM_weights.txt",
                    row.names = F,col.names = F,sep=",")
        
        # gene_comparisons(rval_df)
      }
      # Table
      sdf <- rval_df %>% group_by(community_name) %>% summarize(m = mean(rvals_abs)) %>% arrange(-m)
      sdf_signed <- rval_df %>% group_by(community_name) %>% summarize(m = mean(rvals)) %>% arrange(-m)
      
      # Summary bar plot
      if (plots == T) {
        print(knitr::kable(sdf, caption = "Mean nodal sum(W) per community", floating.environment="sidewaystable"))
        community_plot <- ggplot(data = sdf, aes(x = reorder(community_name,-m), y = m, color = community_name,fill = community_name)) + 
          geom_col(show.legend = F) + scale_fill_brewer(type = "qual",palette = "Dark2") + scale_color_brewer(type = "qual",palette = "Dark2") +
          ylab("Mean weight") + xlab('community') + theme(axis.text.x = element_text(angle = 90,hjust = 1))
        print(community_plot)
        
        print(knitr::kable(sdf_signed, caption = "Mean nodal signed sum(W) per community", floating.environment="sidewaystable"))
        community_plot_signed <- ggplot(data = sdf_signed, aes(x = reorder(community_name,-m), y = m, color = community_name,fill = community_name)) + 
          geom_col(show.legend = F) + scale_fill_brewer(type = "qual",palette = "Dark2") + scale_color_brewer(type = "qual",palette = "Dark2") +
          ylab("Mean weight") + xlab('community')+ theme(axis.text.x = element_text(angle = 90,hjust = 1))
        print(community_plot_signed)
        
        p <- ggplot(data = rval_df,
                    aes(y=rvals_abs,
                        x=ordered(names,levels=rev(names)),
                        fill = community_name,
                        color = community_name)
        ) + 
          geom_bar(stat="identity",show.legend = F) + ylab("W (node sum)") + 
          scale_fill_brewer(type = "qual",palette = "Dark2") + scale_color_brewer(type = "qual",palette = "Dark2") +
          theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm"))+ 
          coord_flip() 
        
        # matrix plot
        g <- expand.grid(x=rownames(commCC),y=colnames(commCC))
        g$vals = c(commCC)
        g$x<-ordered(g$x,levels=colnames(commCC))
        g$y<-ordered(g$y,levels=colnames(commCC))
        grev<-g%>%map_df(rev)
        # Weight matrix
        # matplot <- ggplot(data=grev,
        #                   aes(x=x,
        #                       y=ordered(y,levels=rev(levels(y))),
        #                       fill=vals))+ 
        #   geom_tile(show.legend = T) +
        #   scale_fill_gradient2(low = "blue",high = "red",mid="white",midpoint = 0) + labs(fill = "W (abs val)") +
        #   theme(axis.text.x = element_text(angle = 90,hjust = 0),axis.title = element_blank(),legend.position = "bottom",axis.text = element_text(size = 5),plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        #   scale_x_discrete(position = "top")  + ggtitle(atlasname)
        
        # comboplot_top<-cowplot::plot_grid(plotlist = list(matplot,p),align = "hv",axis = "bt",rel_widths = c(1,.25))
        # comboplot_bottom <- cowplot::plot_grid(plotlist = list(chord_plot,community_plot),align = "hv",axis = "bt",rel_widths = c(1,1))
        # comboplot<-cowplot::plot_grid(comboplot_top,community_plot,ncol = 1,rel_heights = c(1,.33))
        # cowplot::save_plot(sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_matrix.png",atlasname), 
                           # comboplot, base_height = 14, base_width = 10)
        
        if (GABA == T) {
          cat("GABRA gene analysis....\n")
          corr_vec <- vector(mode="double",length = 6)
          for (n in c(1,2,3,4,5,6)) {
            GABA_vals <- read.table(sprintf(
              '/cbica/projects/alpraz_EI/input/atlases/GABRA%d/GABRA%d_%s_threshold_0.95_vals.csv',n,n,atlasname),
              header = T)
            
            GABA_mat_labels <- rval_df %>% left_join(GABA_vals,by = c("nums"="label"))
            GABA_mat_labels$node_names <- GABA_mat_labels$names
            GABA_mat_labels$GABRA <- GABA_mat_labels$value
            colnames(GABA_mat_labels)<- gsub(colnames(GABA_mat_labels),pattern="GABRA",replacement = sprintf("GABRA%d",n))
            GABA_mat_labels$SVM_weights <- GABA_mat_labels$rvals
            GABA_mat_labels$absolute_SVM_weights <- GABA_mat_labels$rvals_abs
            brainsmash(GABA_mat_labels,c(sprintf("GABRA%d",n),"SVM_weights"))
            brainsmash(GABA_mat_labels,c(sprintf("GABRA%d",n),"absolute_SVM_weights"))
            # gaba_plot <- region_plots(df = GABA_mat_labels,var_name = "rvals",xlabel = sprintf("GABRA%d",n))
            # print(gaba_plot)
            cortex = GABA_mat_labels %>% filter(!(community_name%in%c("AAL_subcortex")))
            corr_vec[n] <- cor(cortex$SVM_weights,cortex$GABRA)
          }
          if (file.exists('surrogate_brainmap_corrs_GABRA3_map.txt')) {
            x=data.table::rbindlist(lapply(list.files(pattern = glob2rx("surr*brainmap*GAB*.txt")), data.table::fread),idcol="GABRA")
            names(x) <- c("GABRA", "r")
            x<-x %>%
              mutate(GABRA = factor(GABRA,levels = c(1,2,3,4,5,6),labels = c("GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6")),
                     R2=r^2,
                     source = "BrainSMASH null")
            corr_df <- data.frame(GABRA = factor(c("GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6")),
                                  r = corr_vec,
                                  R2 = corr_vec^2,
                                  source = "Observed")
            x <- rbind(x,corr_df)

            rc_plot <- ggplot(data=x,aes(x = GABRA,y = R2,color = source)) + 
              geom_violin(show.legend = TRUE,alpha = .9,fill="black") +
              # geom_point(aes(color=source,size=source,alpha=source),position = position_jitter(width = .1),show.legend = TRUE) +
              # scale_color_manual(values = c("black","red"),labels = c("BrainSMASH null","Observed")) +
              # scale_size_manual(values = c(.25,3))+
              # scale_alpha_manual(values = c(.25,1))+
              ylab(bquote("Spatial relationship" ~ (italic(R)^2)))+theme(legend.title = element_blank(),axis.title.x = element_blank())
            
            rc_plot <- rc_plot + geom_point(data=corr_df,aes(x=as.factor(GABRA),y=R2,color="red"),size=3,show.legend = TRUE) +
              scale_color_manual(values = c("black","red"),labels = c("BrainSMASH null","Observed")) +
              theme(legend.title = NULL)
            print(rc_plot)
            
            # rc_plot <- ggplot(data=x,aes(y = GABRA,x = R2,color = source)) + 
            #   geom_density_ridges(show.legend = TRUE) +
            #   # geom_point(aes(color=source,size=source,alpha=source),position = position_jitter(width = .1),show.legend = TRUE) +
            #   # scale_color_manual(values = c("black","red"),labels = c("BrainSMASH null","Observed")) +
            #   # scale_size_manual(values = c(.25,3))+
            #   # scale_alpha_manual(values = c(.25,1))+
            #   ylab(bquote("Spatial relationship" ~ (italic(R)^2)))+theme(legend.title = element_blank(),axis.title.x = element_blank())
            # 
            # rc_plot <- rc_plot + geom_point(data=corr_df,aes(y=as.factor(GABRA),x=R2,color="red"),size=3,show.legend = TRUE) +
            #   scale_color_manual(values = c("black","red"),labels = c("BrainSMASH null","Observed")) +
            #   theme(legend.title = NULL)
            # 
          }
          cat("GABRA gene analysis complete \n")
        }
        
        if (transmodality == T){
          if (atlasname == "schaefer200x7_aal") {
            trans_labels <- read.csv('/cbica/projects/alpraz_EI/input/atlases/Schaef200_transmodality7.csv',header = F)
          } else{
            trans_labels <- read.csv('/cbica/projects/alpraz_EI/input/atlases/Schaef400_transmodality7.csv',header = F)
          }
          
          colnames(trans_labels) <- c('name','transmodality')
          trans_labels$name <- gsub(trans_labels$name,pattern = "7Networks_",replacement = "")
          trans_labels$nums = as.numeric(rownames(trans_labels))

          mat_labels <- mat_labels %>%
            left_join(trans_labels,by = "nums")
          
          rval_df_trans <- rval_df %>% select(-node_community,-nums,-community_name)%>%
            left_join(mat_labels,by = c("names" = "node_names"),keep=TRUE)

          #Brainsmash
          rval_df_trans$SVM_weights <- rval_df_trans$rvals
          rval_df_trans$absolute_SVM_weights <- rval_df_trans$rvals_abs
          brainsmash(rval_df_trans,c("transmodality","SVM_weights"))
          brainsmash(rval_df_trans,c("transmodality","absolute_SVM_weights"))
          
          colnames(newCC)=mat_labels$node_names
          rownames(newCC)=mat_labels$node_names
          # Arrange by community membership
          orderCC <- newCC[order(mat_labels$transmodality),order(mat_labels$transmodality)]
          newCC <- orderCC
          
          # matrix plot
          g <- expand.grid(x=rownames(newCC),y=colnames(newCC))
          g$vals = c(newCC)
          g$x<-ordered(g$x,levels=colnames(commCC))
          g$y<-ordered(g$y,levels=colnames(commCC))
          grev<-g%>%map_df(rev)
          # Weight matrix
          matplot <- ggplot(data=grev,
                            aes(x=x,
                                y=ordered(y,levels=rev(levels(y))),
                                fill=vals))+ 
            geom_tile(show.legend = T) +
            scale_fill_gradient2(mid = "white",high = "red",low="blue",midpoint = 0) + labs(fill = "W (abs val)") +
            theme(axis.text.x = element_text(angle = 90,hjust = 0),axis.title = element_blank(),legend.position = "bottom",axis.text = element_text(size = 5),plot.margin = unit(c(0, 0, 0, 0), "cm")) +
            scale_x_discrete(position = "top")  + ggtitle(paste0('transmodality ',atlasname))
          
          r = cor.test(rval_df_trans$transmodality,rval_df_trans$rvals)
          dotplot <- ggplot(data=rval_df_trans,
                            aes(x=transmodality,
                                y=rvals,
                                color = community_name)) +
            geom_point()+
            geom_smooth(method = "lm",aes(x=transmodality,y=rvals,color=NULL)) +
            annotate(x=mean(rval_df_trans$transmodality,na.rm=T),y=mean(rval_df_trans$rvals),geom = "text",label=sprintf("r = %1.2f,p = %1.3f",r$estimate,r$p.value))+
            ggtitle("Region sum W by transmodality")+ylab("Region global SVM weight")
          tlegend <- get_legend(dotplot +theme(legend.position = "bottom",legend.title = element_blank(),legend.background = element_blank()))
          dotplot <- dotplot + theme(legend.position = "none")
          # print(dotplot)
          #absolute value
          r = cor.test(rval_df_trans$transmodality,rval_df_trans$rvals_abs)
          dotplot_abs <- ggplot(data=rval_df_trans,
                                aes(x=transmodality,
                                    y=rvals_abs,
                                    color = community_name)) +
            geom_point()+
            geom_smooth(method = "lm",aes(x=transmodality,y=rvals_abs,color=NULL)) +
            annotate(x=mean(rval_df_trans$transmodality,na.rm=T),y=mean(rval_df_trans$rvals_abs),geom = "text",label=sprintf("r = %1.2f,p = %1.3f",r$estimate,r$p.value))+
            ggtitle("Region sum absolute W by transmodality")+theme(legend.position = "none")+ylab("Region global SVM Weight")
          dotplots_top <- cowplot::plot_grid(dotplot,dotplot_abs,nrow = 1,align = "h",axis = "lr")
          dotplots <- cowplot::plot_grid(dotplots_top,tlegend,nrow = 2,align = "h",rel_heights = c(1,.3))
          print(dotplots)
          
          g <- g %>% left_join(mat_labels%>% select(node_names,transmodality),by = c("x"="node_names"))
          g <- g %>% left_join(mat_labels%>% select(node_names,transmodality),by = c("y"="node_names"),suffix = c("x","y"))

          comboplot<-cowplot::plot_grid(plotlist = list(matplot,dotplot),ncol=1,align = "hv",axis = "lr",rel_heights = c(1,.66))
          ggdraw(comboplot)
          # print(comboplot)
          cowplot::save_plot(sprintf("/cbica/projects/alpraz_EI/output/drug_classification/transmodality_%s_matrix.png",atlasname), 
                             comboplot, base_height = 12, base_width = 8)
        }
        
        return_obj <- list(feat_mat = newCC,community_ranks = sdf)
      } else {
        if (returnRankingsOnly == T) {
          return_obj <- sdf
        } else {
          return_obj <- list(feat_mat = newCC,community_ranks = sdf)
        }
      }
    } else{
      return_obj <- newCC
    }
    
  } else {
    mat_labels <- data.frame(nums = colnames(templateCC))
    mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
    mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
    colnames(newCC)=mat_labels$V2
    rownames(newCC)=mat_labels$V2
    return_obj <- list(feat_mat = newCC)
  }
  
  return(return_obj)
}

mat2brain <- function(CC,atlasname,filename){
  
  #load atlas
  fname <- sprintf("/cbica/projects/alpraz_EI/input/atlases/%s_threshold_0.95.nii.gz",atlasname)
  atlas <- readnii(fname)
  atlas_vec <- c(atlas)
  
  # Get labels
  fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK_GSR/13783/1986/%s_CC_GSR_000.netcc",atlasname)
  template_CC <- as.matrix(read.table(fname,header = T,skip = 4))
  labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_labels.csv',atlasname),header=F)
  if (atlasname %in% c("schaefer200x7_aal","schaefer400x7_aal","gordon333_aal")) {
    community_labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_node_labels.csv',
                                         gsub(x = atlasname,pattern = "_aal",replacement = ""))
                                 ,header=T)
    

    mat_labels <- data.frame(nums = colnames(template_CC))
    mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
    mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
    mat_labels <- mat_labels %>% left_join(community_labels,by = c("nums"="index"))
    mat_labels <- mat_labels %>%
      mutate(community_name = as.character(community_name),V2 = as.character(V2))
    # mat_labels$community_name[is.na(mat_labels$community_name)]=mat_labels$V2[is.na(mat_labels$community_name)]
    mat_labels$community_name[is.na(mat_labels$community_name)]="AAL_subcortex"
    
    # Calculate values to map on brain
    vals <-rowSums(CC)
    vals <- data.frame(vals = vals,names = names(vals))
    ## Get key column
    idx <- sapply(mat_labels, function(x) sum(x%in%vals$names)/length(vals$names))
    matching_col<-colnames(mat_labels)[idx==T][1] # get first match if multiple
    mat_labels$names <- mat_labels[,matching_col]
    label_df <- mat_labels %>% 
      left_join(vals, by = "names")
    
    # Map values on brain
    voxel_df <- as.data.frame(atlas_vec) %>% left_join(label_df,by = c("atlas_vec"="nums"))
    label_vec <- as.matrix(voxel_df$vals)
    label_matrix <- array(label_vec, dim = dim(atlas))
    nifti_image <- niftiarr(atlas,label_matrix)
    cat(sprintf("\n Writing /cbica/projects/alpraz_EI/output/drug_classification/%s.nii.gz ......\n",filename))
    writenii(nifti_image,filename = sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s.nii.gz",filename))
    cat("done\n")
    
  } else {
    mat_labels <- data.frame(nums = colnames(template_CC))
    mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
    mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
    colnames(newCC)=mat_labels$V2
    rownames(newCC)=mat_labels$V2
  }
  
}

vec2surface <- function(schaefer_df,val_name,file_prefix){
  cat("\nWriting Surface Files..\b")
  # Left side
  left_idx <- read.csv('/cbica/projects/alpraz_EI/input/atlases/400L_idx.csv',header = F)
  colnames(left_idx) <- "IDX"
  left_table <- read.csv('/cbica/projects/alpraz_EI/input/atlases/400L.csv',header = T)
  colnames(left_table) <- c("Var1","Var2","Var3","Var4","IDX","node_names")
  #Tidy up some of the names
  left_table$node_names <- gsub(left_table$node_names,pattern = "7Networks_",replacement = "")
  left_table$node_names= gsub(left_table$node_names,pattern = "Default_pCun",replacement = "Default_")
  left_table$node_names= gsub(left_table$node_names,pattern = "FrOperIns",replacement = "FrOper")

  schaefer_df$value <- as.numeric(schaefer_df[,val_name])

  schaefer_df_left <- left_table %>%
    left_join(schaefer_df,by = c("node_names"="names")) %>%
    select(node_names,IDX,value)
  
  left_idx_values <- left_idx %>%
    left_join(schaefer_df_left,by="IDX")
  # left_idx_values$value[is.na(left_idx_values$value)]=0
  left_output <- left_idx_values$value
  write.table(left_output,sprintf("/cbica/projects/alpraz_EI/output/surface_files/%s_LH.csv",file_prefix),row.names = F,col.names = F)

  # Right side
  right_idx <- read.csv('/cbica/projects/alpraz_EI/input/atlases/400R_idx.csv',header = F)
  colnames(right_idx) <- "IDX"
  right_table <- read.csv('/cbica/projects/alpraz_EI/input/atlases/400R.csv',header = T)
  colnames(right_table) <- c("Var1","Var2","Var3","Var4","IDX","node_names")
  #Tidy up some of the names
  right_table$node_names <- gsub(right_table$node_names,pattern = "7Networks_",replacement = "")
  right_table$node_names= gsub(right_table$node_names,pattern = "Default_pCun",replacement = "Default_")
  right_table$node_names= gsub(right_table$node_names,pattern = "FrOperIns",replacement = "FrOper")
  right_table$node_names= gsub(right_table$node_names,pattern = "Default_PFCdPFCm",replacement = "Default_PFCm")
  
  schaefer_df$value <- as.numeric(schaefer_df[,val_name])
  schaefer_df_right <- right_table %>%
    right_join(schaefer_df,by = c("node_names"="names")) %>%
    select(node_names,IDX,value)

  right_idx_values <- right_idx %>%
    left_join(schaefer_df_right,by="IDX")
  # right_idx_values$value[is.na(right_idx_values$value)]=0
  right_output <- right_idx_values$value
  write.table(right_output,sprintf("/cbica/projects/alpraz_EI/output/surface_files/%s_RH.csv",file_prefix),row.names = F,col.names = F)
  cat("done\n")
}


region_plots <- function(df,var_name,xlabel){
  ## signed plots ##
  cat(sprintf('\n## Signed W for %s\n',xlabel))
  # all regions
  r=cor.test(df[,var_name],df$value,method = "pearson")
  
  all_plot <- ggplot(data = df,aes_string(x="value",y=var_name,color="community_name")) + 
    geom_point() +geom_smooth(method = "lm") + geom_smooth(aes_string(x="value",y=var_name),color="black",method = "lm") +
    xlab(xlabel) + ylab("mean W (signed)")+
    annotate(y=mean(df[,var_name]),x=mean(df$value),geom = "text",label=sprintf("r = %1.2f, p = %1.3f",r$estimate,r$p.value))+
    ggtitle("signed W")
  legend_all = get_legend(all_plot+theme(legend.position = "right",legend.background=element_blank(),legend.box.margin = margin(0,0,0,12),legend.title = element_blank()))
  all_plot <- all_plot+theme(legend.position = "none")
  
  # cortex
  cortex = df %>% filter(!(community_name%in%c("AAL_subcortex")))
  r=cor.test(cortex[,var_name],cortex$value,method = "pearson")
  cortical_plot <- ggplot(data = cortex,aes_string(x="value",y=var_name,color="community_name")) +
    geom_point() +xlab(xlabel) +
    ylab("mean W (signed)")+geom_smooth(method = "lm")+geom_smooth(aes_string(x="value",y=var_name),color="black",method = "lm") +
    annotate(y=mean(cortex[,var_name]),x=mean(cortex$value),geom = "text",label=sprintf("r = %1.2f, p = %1.3f",r$estimate,r$p.value))+
    ggtitle("cortical SVM W")
  legend_cortex = get_legend(cortical_plot+theme(legend.position = "right",legend.background=element_blank(),legend.box.margin = margin(0,0,0,12),legend.title = element_blank()))
  cortical_plot <- cortical_plot+theme(legend.position = "none")
  
  #subcortex
  subcortex = df %>% filter(community_name%in%c("AAL_subcortex"))
  r=cor.test(subcortex[,var_name],subcortex$value,method = "pearson")
  subcortical_plot <- ggplot(data = subcortex,aes_string(x="value",y=var_name,color="community_name")) +
    geom_point() +xlab(xlabel) +
    ylab("mean W (signed)")+geom_smooth(method = "lm")+
    annotate(y=mean(subcortex[,var_name]),x=mean(subcortex$value),geom = "text",label=sprintf("r = %1.2f, p = %1.3f",r$estimate,r$p.value))+
    ggtitle("sucortical SVM W")+theme(legend.position = "right",legend.background=element_blank())
  legend_subcortex = get_legend(subcortical_plot+theme(legend.position = "right",legend.background=element_blank(),legend.box.margin = margin(0,0,0,12),legend.title = element_blank()))
  subcortical_plot <- subcortical_plot+theme(legend.position = "none")
  
  ## Absolute value ##
  cat(sprintf('\n## Absolute W for %s\n',xlabel))
  # all regions
  var_name = sprintf("%s_abs",var_name)
  r=cor.test(df[,var_name],df$value,method = "pearson")
  abs_all_plot <- ggplot(data = df,aes_string(x="value",y=var_name,color="community_name")) + 
    geom_point() +geom_smooth(method = "lm") + geom_smooth(aes_string(x="value",y=var_name),color="black",method = "lm") +
    xlab(xlabel) + ylab("mean abs(W)")+
    annotate(y=mean(df[,var_name]),x=mean(df$value),geom = "text",label=sprintf("r = %1.2f, p = %1.3f",r$estimate,r$p.value))+
    ggtitle("absolute W")+theme(legend.position = "none")
  
  # cortex
  r=cor.test(cortex[,var_name],cortex$value,method = "pearson")
  abs_cortical_plot <- ggplot(data = cortex,aes_string(x="value",y=var_name,color="community_name")) +
    geom_point() +xlab(xlabel) +
    ylab("mean abs(W)")+geom_smooth(method = "lm")+geom_smooth(aes_string(x="value",y=var_name),color="black",method = "lm") +
    annotate(y=mean(cortex[,var_name]),x=mean(cortex$value),geom = "text",label=sprintf("r = %1.2f, p = %1.3f",r$estimate,r$p.value))+
    ggtitle("cortical absolute W")+theme(legend.position = "none")
  
  #subcortex
  r=cor.test(subcortex[,var_name],subcortex$value,method = "pearson")
  abs_subcortical_plot <- ggplot(data = subcortex,aes_string(x="value",y=var_name,color="community_name")) +
    geom_point() +xlab(xlabel) +
    ylab("mean abs(W)")+geom_smooth(method = "lm")+
    annotate(y=mean(subcortex[,var_name]),x=mean(subcortex$value),geom = "text",label=sprintf("r = %1.2f, p = %1.3f",r$estimate,r$p.value))+
    ggtitle("sucortical absolute W")+theme(legend.position = "none")
  
  top_row <- cowplot::plot_grid(all_plot,abs_all_plot,legend_all,nrow=1,rel_widths = c(1,1,.333))
  cortex_row <- cowplot::plot_grid(cortical_plot,abs_cortical_plot,legend_cortex,nrow=1,rel_widths = c(1,1,.333))
  subcortex_row <- cowplot::plot_grid(subcortical_plot,abs_subcortical_plot,legend_subcortex,nrow=1,rel_widths = c(1,1,.333))
  combined_plot <- cowplot::plot_grid(top_row,cortex_row,subcortex_row,ncol = 1)
  return(combined_plot)
}

write_brainsmash_files <- function(df,distance_matrix,idx,fname,var_names){
  this_distance_matrix <- distance_matrix[idx,idx]
  write.table(this_distance_matrix,file = sprintf("/cbica/projects/alpraz_EI/output/brainsmash/input_files/schaefer400x7_aal_threshold_0.95_%s_distMat.csv",fname),
              sep = ",",row.names = F,col.names = F)
  brain_map1 <- df[idx,var_names[1]]
  write.table(brain_map1,file = sprintf("/cbica/projects/alpraz_EI/output/brainsmash/input_files/schaefer400x7_aal_threshold_0.95_%s_%s.csv",fname,var_names[1]),
              sep = ",",row.names = F,col.names = F)
  brain_map2 <- df[idx,var_names[2]]
  write.table(brain_map2,file = sprintf("/cbica/projects/alpraz_EI/output/brainsmash/input_files/schaefer400x7_aal_threshold_0.95_%s_%s.csv",fname,var_names[2]),
              sep = ",",row.names = F,col.names = F)
  
}

brainsmash <- function(df,var_names){
  distance_matrix <- read.table("/cbica/projects/alpraz_EI/input/atlases/schaefer400x7_aal_threshold_0.95_distMat.txt")
  distance_matrix_expanded <- squareform(distance_matrix$V1)
  
  distance_matrix_labels <- read.table("/cbica/projects/alpraz_EI/input/atlases/schaefer400x7_aal_threshold_0.95_roi_order.csv",
                                       header = T,sep = ",",row.names = 1)
  colnames(distance_matrix_labels) <- "labels"
  labeled_brain_map <- left_join(distance_matrix_labels,df,by = c("labels"="nums"),keep=TRUE) 

  # Left cortex
  left_parcel_numbers <- labeled_brain_map$nums[grepl(x=labeled_brain_map$node_names,pattern = "LH")]
  left_parcel_idx <- match(left_parcel_numbers,labeled_brain_map$labels)#this step is probably not necessary but is for abundance of caution for indexing.
  
  write_brainsmash_files(df = labeled_brain_map,distance_matrix = distance_matrix_expanded,idx = left_parcel_idx,fname = "left",var_names = var_names)
  
  # Right cortex
  right_parcel_numbers <- labeled_brain_map$nums[grepl(x=labeled_brain_map$node_names,pattern = "RH")]
  right_parcel_idx <- match(right_parcel_numbers,labeled_brain_map$labels)#this step is probably not necessary but is for abundance of caution for indexing.
  
  write_brainsmash_files(df = labeled_brain_map,distance_matrix = distance_matrix_expanded,idx = right_parcel_idx,fname = "right",var_names = var_names)
  
  # All cortex
  cortex_parcel_numbers <- labeled_brain_map$nums[!is.na(labeled_brain_map$node_community)]
  cortex_parcel_idx <- match(cortex_parcel_numbers,labeled_brain_map$labels)#this step is probably not necessary but is for abundance of caution for indexing.
  
  write_brainsmash_files(df = labeled_brain_map,distance_matrix = distance_matrix_expanded,idx = cortex_parcel_idx,fname = "cortex",var_names = var_names)
  
  # Subcortex
  subcortex_parcel_numbers <- labeled_brain_map$nums[labeled_brain_map$community_name=="AAL_subcortex"]
  subcortex_parcel_idx <- match(subcortex_parcel_numbers,labeled_brain_map$labels)#this step is probably not necessary but is for abundance of caution for indexing.
  
  write_brainsmash_files(df = labeled_brain_map,distance_matrix = distance_matrix_expanded,idx = subcortex_parcel_idx,fname = "subcortex",var_names = var_names)

}

gene_comparisons <- function(df){
  gene_labels <- read.table('/cbica/projects/alpraz_EI/input/atlases/genes.csv',sep = ",",header = T)%>%
    mutate(entrezgene_id = as.character(entrezgene_id))
  load('/cbica/projects/alpraz_EI/input/atlases/AHBA_vertex_Schaefer4007order.RData')
  # rh_genes=rh_genes_mirr
  # lh_genes = lh_genes_mirr
  rm(rh_genes_mirr,lh_genes_mirr)
  #Tidy up left hemi names
  colnames(lh_genes)<- gsub(colnames(lh_genes),pattern = "7Networks_",replacement = "")
  colnames(lh_genes)<- gsub(colnames(lh_genes),pattern = "Default_pCun",replacement = "Default_")
  colnames(lh_genes)<- gsub(colnames(lh_genes),pattern = "FrOperIns",replacement = "FrOper")
  
  #Tidy up right hemi names
  colnames(rh_genes)<- gsub(colnames(rh_genes),pattern = "7Networks_",replacement = "")
  colnames(rh_genes)<- gsub(colnames(rh_genes),pattern = "Default_pCun",replacement = "Default_")
  colnames(rh_genes)<- gsub(colnames(rh_genes),pattern = "FrOperIns",replacement = "FrOper")
  colnames(rh_genes)<- gsub(colnames(rh_genes),pattern = "Default_PFCdPFCm",replacement = "Default_PFCm")
  
  genes <- cbind(lh_genes,rh_genes)
  column_idx <- colnames(genes) %in% df$names
  genes_idxed <- genes[,column_idx]
  genes_t <- as.data.frame(t(genes_idxed))
  genes_t$names <- rownames(genes_t)
  
  genes_t_joined <- genes_t %>% left_join(df %>%
                                            select(names,rvals,rvals_abs),
                                          by = "names")
  rm(lh_genes,rh_genes,genes_idxed,genes_t)
  
  genes <- as.matrix(genes_t_joined %>%
                       select(-rvals,-rvals_abs,-names))
  W_vals <- as.matrix(genes_t_joined %>%
                        select(rvals,rvals_abs))
  cc <- cor(genes,W_vals[,"rvals"])
  cc_abs <- cor(genes,W_vals[,"rvals_abs"])
  cc_df <- data.frame(cc=cc,label=row.names(cc)) %>% 
    arrange(cc) %>%
    left_join(gene_labels,by = c("label"="entrezgene_id"))%>%
    mutate(hgnc_symbol=as.character(hgnc_symbol))
  cc_abs_df <- data.frame(cc=cc_abs,label=row.names(cc)) %>% 
    arrange(cc)%>%
    left_join(gene_labels,by = c("label"="entrezgene_id"))
  GABRA_corrs <- cc_df$cc[str_detect(cc_df$hgnc_symbol,"GABRA")]
  GABRA_corrs <- GABRA_corrs[!is.na(GABRA_corrs)]
  GABRA_labs <- cc_df$hgnc_symbol[str_detect(cc_df$hgnc_symbol,"GABRA")]
  GABRA_labs <- as.character(GABRA_labs[!is.na(GABRA_labs)])
  GABRA_ps <- sapply(GABRA_corrs,FUN = function(x) sum(abs(cc_df$cc)>abs(x))/length(cc_df$cc))
  print(rbind(GABRA_labs,GABRA_corrs,GABRA_ps))

  outplot<-ggplot(cc_df,aes(x=cc)) + geom_histogram() + geom_vline(xintercept = GABRA_corrs) +
    annotate(geom="text",x=GABRA_corrs,y=10,label = GABRA_labs)
  
  print(outplot)
  print(knitr::kable(cc_df%>%head(50)))
  print(knitr::kable(cc_df%>%tail(50)))

}

ROC_curve <- function(DecisionValues, labels){
  # Decision values is nx1
  # Labels is nx1
  
  # N.B.
  # Drug (class 0) is assigned as +1 by libsvm and placebo (class 1) is assigned as 0 by libsvm default. 
  # Adjust the true labels and predicted labels to match that here so that the decision values make sense.
  labels <- labels*-1+1
  
  P <- sum(labels == 1)
  N <- sum(labels == 0)
  Sorted_DecisionValues <- sort(unique(DecisionValues), decreasing = FALSE)
  numDecisionValues <- length(Sorted_DecisionValues)
  
  TP_Array <- vector(mode = "numeric",length = numDecisionValues)
  FP_Array <- vector(mode = "numeric",length = numDecisionValues)
  Accuracy_Array = vector(mode = "numeric",length = numDecisionValues)
  for (i in 1:numDecisionValues){
    thisCutoff <- Sorted_DecisionValues[i]
    thisPredictedLabels <- as.numeric(DecisionValues>thisCutoff)
    detections <- thisPredictedLabels==1

    TP <- sum(labels[detections] == thisPredictedLabels[detections])
    TPR <- TP/P
    FP <- sum(labels[detections]!=thisPredictedLabels[detections])
    FPR <- FP/N
    
    TP_Array[i] <- TPR
    FP_Array[i] <- FPR
    
    Accuracy_Array[i] = (TP + N - FP) / (P + N)
  }
  
  LargestAccuracy = max(Accuracy_Array)
  ROC_output <- data.frame(TPR =TP_Array,FPR=FP_Array,Accuracy = Accuracy_Array)
  ROC_output <- ROC_output%>%arrange(TPR,FPR)
  
  #AUC
  dFPR <- c(0,diff(ROC_output$FPR))
  dTPR <- c(0,diff(ROC_output$TPR))
  AUC <- sum(ROC_output$TPR * dFPR) + sum(dTPR * dFPR)/2
  
  # Plot 
  MaxAccuracyIdx <- which.max(ROC_output$Accuracy)
  roc_plot <- ggplot(ROC_output,aes(x=FPR,y=TPR))+geom_point() +
    geom_line() +
    geom_point(x=ROC_output$FPR[MaxAccuracyIdx],y=ROC_output$TPR[MaxAccuracyIdx],size=3,aes(color="red"))+
    scale_color_manual(values="red",labels=sprintf("Max Accuracy = %1.3f",LargestAccuracy))+
    theme(legend.title = element_blank(),legend.position = c(.5,.25),legend.background = element_blank(),legend.justification = c("left","top"))+
    annotate("text",x=.58,y=.2,label = sprintf("AUC = %1.3f",AUC),hjust=0,size= theme_get()$text[["size"]]/4)+
    xlab('False positive rate (1-specificity)')+ylab("True positive rate (sensitivity)")
  print(roc_plot)
}

display_results <- function(atlasname,GSR="GSR",classifier="svm",perm_results=F,result_fname=NULL,results=NULL){
  cat(sprintf("\n## Displaying classification results for %s atlas (classifier = %s)\n",atlasname,classifier))
  
  # df_acc <- read.csv(sprintf('/cbica/projects/alpraz_EI/output/drug_classification/%s_results_%s.csv',GSR,classifier))
  # acc_plot <- ggplot(data = df_acc,aes(x=num_features,y=accuracy,color = atlas,group=atlas))+
  #   geom_point()+
  #   ggtitle(sprintf("Accuracy for %s %s classification",GSR,classifier))
  # print(acc_plot)
  if (is.null(results)) {
    if (!is.null(result_fname)) {
      results <- readRDS(result_fname)
    } else {
      if (perm_results == T){
        results <- readRDS(sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_1_permute_on_results.rds",
                                   atlasname,GSR,classifier))
      }else {
        results <- readRDS(sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_all_%s_1_permute_off_results.rds",
                                   atlasname,GSR,classifier))
      }
    }
  }
  
  b <- results[[1]]
  cat("Results for exact binomial test:\n")
  cat(sprintf('accuracy = %1.3f',b$accuracy))
  print(b)
  
  # Look at W coefs
  W <- results[[2]]
  if (atlasname=="schaefer400x7_aal") {
    feat_mat_obj <- feat2mat(W,atlasname,community_summary = T,GABA=T,transmodality = T,surface = T)
  } else{
    feat_mat_obj <- feat2mat(W,atlasname,community_summary = T,GABA=T,transmodality = T)
  }
  Wmap <- feat_mat_obj$feat_mat
  comm_ranks <- feat_mat_obj$community_ranks
  # c1 <- corrplot::corrplot(Wmap,is.corr = F, method = "color",na.label = "square",title = sprintf("%s W coefficients",atlasname),
  #                            diag = F,tl.cex = .2)
  # knitr::include_graphics(sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_matrix.png",atlasname))
  # mat2brain(CC = abs(Wmap),atlasname = atlasname, filename = sprintf("%s_%s_%s",atlasname,GSR,classifier))

  if ("perm_W" %in% names(b)) {
    # show permuted W significance
    b <- results[[1]]
    cat("Results for alpraz (exact binomial test):\n")
    cat(sprintf("\n Number of features: %d\n",results[[3]]))
    print(b$estimate)
    cat(sprintf("p = %1.3f\n",b$p.value))
    cat(sprintf("Permutation p = %1.3f\n",b$perm_p))
    perm_plot <- ggplot(data = data.frame(perm_acc=t(b$perm_results[[1]])),aes(x = perm_acc)) +
      geom_histogram()+geom_vline(xintercept = b$accuracy) + 
      geom_label(x = b$estimate,y = 100,label=paste0("Observed\n p = ",as.character(round(b$perm_p,digits = 4))))+
      xlab("Classification Accuracy")+ylab("Number of Draws")+
      ggtitle("Permutation Test")+coord_cartesian(clip = "off")
      
    print(perm_plot)
    
    pred_data <- b$pred_data
    ROC_curve(pred_data$decisionValues,pred_data$drug)
    
    
    # permW <- b$perm_W
    # if (file.exists(sprintf("%s_permutation_communities.rds",atlasname))){
    #   perm_comm_summary <- readRDS('schaefer200x7_aal_permutation_communities.rds')
    # } else{
    #   perm_comm_summary <- perm_communities(permW,atlasname)
    #   saveRDS(perm_comm_summary,file = sprintf("%s_permutation_communities.rds",atlasname))
    # }
    # perm_comm_summary<- lapply(perm_comm_summary,FUN = function(x){arrange(x,community_name)})
    # summary_df <- as.data.frame(sapply(perm_comm_summary, `[`, i =, j = 2))
    # summary_df <- summary_df[match(comm_ranks$community_name,perm_comm_summary[[1]]$community_name),]
    # perm_comm_summary = perm_comm_summary[match(comm_ranks$community_name,perm_comm_summary$community_name),]
    # c_test <- sapply(perm_comm_summary, FUN = function(x) {
    #   x>comm_ranks$m
    # })
    # rowMeans(as.data.frame(c_test) %>% select(contains("W")))
    # rowSums(as.data.frame(c_test) %>% select(contains("W")))
    # long_summary <- as.data.frame(perm_comm_summary) %>% pivot_longer(cols = select(contains("W")))
    
    # # Feature level sig test
    # w_sig <- b$perm_W_sig[[1]]
    # w_sig_thresh = w_sig<.05
    # feat_mat_obj<-feat2mat(w_sig_thresh,atlasname = atlasname,community_summary = T)
    # wsig_mat <- feat_mat_obj$feat_mat
  }
  else{
    pred_data <- b$pred_data[[1]] #Grabbing first one just to get the drug labels, subid, sesid.
    dec_folds <- data.table::rbindlist(b$pred_data,idcol = "fold")%>%
      pivot_wider(names_from = "fold","values_from"="decisionValues",id_cols = c("subid","sesid")) %>% 
      group_by(subid,sesid)
    pred_data$mean_dec_vals <- rowMeans(dec_folds[,3:dim(dec_folds)[2]])
    ROC_curve(pred_data$mean_dec_vals,pred_data$drug)
  }
  # return(results)
}


### function to extract derivative, confidence interval, significance, and plot ###
### This works for signle smooth terms and factor-smooth interactions. Does not work for bivariate smooths.
### This part is still under development.
### If you want to plot derivatives for a by-variable factor model, and you want all plots to have the same scaling, see the note below. Right now you will have to manually set the max value (sorry)

get_derivs_and_plot <- function(modobj,smooth_var,low_color=NULL,hi_color=NULL){
  this_font_size = font_size
  if (is.null(low_color)){low_color = "white"}
  if (is.null(hi_color)){hi_color = "grey20"}
  derv<-derivatives(modobj,term=smooth_var)
  derv<- derv %>%
    mutate(sig = !(0 >lower & 0 < upper))
  derv$sig_deriv = derv$derivative*derv$sig
  cat(sprintf("\nSig change: %1.2f - %1.2f\n",min(derv$data[derv$sig==T]),max(derv$data[derv$sig==T])))
  d1<- ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = sig_deriv))
  
  # Set the gradient colors
  if (min(derv$derivative)>0) {
    d1 <- d1 + scale_fill_gradient(low = low_color, high = hi_color,limits = c(0,max(derv$derivative)))
    # If you want all plots to have the same scaling, this code can be used instead-- This is desirable if you have a factor-smooth model.
    ## max_val = .5
    ## scale_fill_gradient(low = low_color,high = hi_color,limits = c(0,max_val),oob = squish)
  } else if (min(derv$derivative)<0 & max(derv$derivative)<0) {
    d1 <- d1 + scale_fill_gradient(low = hi_color, high = low_color,limits = c(min(derv$derivative),0))
  }else {
    d1 <- d1 +
      scale_fill_gradient2(low = "steelblue", midpoint = 0, mid = "white",
                           high = "firebrick",limits = c(min(derv$derivative),max(derv$derivative)))
  }

  d1 <- d1 + 
    labs(x = smooth_var,fill = sprintf("\u0394%s",smooth_var)) + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = this_font_size),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=this_font_size),
          legend.text = element_text(size = this_font_size),
          axis.title = element_text(size = this_font_size),
          legend.key.width = unit(1,"cm"),
          legend.position = "bottom",
          legend.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    guides(fill = guide_colorbar(reverse = F,direction = "horizontal",title.position = "left")) +
    geom_rect(aes(ymin=0,ymax=1,xmin=min(data),xmax=max(data)),color="black",fill="white",alpha = 0)
  return(d1)
}

# Func to visualize model outputs
visualize_model <- function(modobj,smooth_var, int_var = NULL ,group_var = NULL, plabels = NULL,check_diagnostics = F,derivative_plot = F){
  this_font_size = font_size
  if (any(class(modobj)=="gam")) {
    model <- modobj
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  s<-summary(model)
  
  ## Generate custom line plot
  np <- 10000 #number of predicted values
  df = model$model
  
  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  
  if (!is.null(int_var)) {
    # We will produce and interaction plot
    if (!any(grepl(x=as.character(model$formula),pattern = int_var))) {
      warning("int_var not recognized in model formula!")
      return()
    }
    switch (varClasses[int_var],
            "numeric" = {
              q <- quantile(df[,int_var],probs = c(.05,.95)) #pick 10% and 90% to plot
              bigq <- q[[2]]
              smallq <- q[[1]]
              values <- c(bigq,smallq)
              labs <- c(sprintf("high (%1.2f)",bigq),sprintf("low (%1.2f)",smallq))
              
              q <-quantile(rescale(df[,int_var],c(0,1)),probs = c(0,.5,1))
              limit_values <- c(q[[1]],q[[length(q)]])
              midpoint_val <- unname(q[[2]])
              cbar_vals <- unname(q)
              
              theseLabs = rep(values,each = np)
              grad_fill = T
            },
            "factor" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = rep(values,each = np)
              grad_fill = F
            },
            "ordered" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = ordered(rep(values,each = np),levels = values)
              grad_fill = F
            }
    )
    labPred <- data.frame(init = rep(0,np*length(labs)))
    labPred[,int_var] = theseLabs
    labPred$lab = rep(labs,each = np)
    labPred <- labPred[,names(labPred) !="init"]
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else if (v == int_var) {
        next
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    
    thisPred <- thisPred %>% select(-init)
    thisPred <- do.call("rbind", replicate(length(labs), thisPred, simplify = FALSE))
    
    pred <- cbind(labPred,thisPred)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    if (!is.null(group_var)) {
      pred[,group_var] = NA #these columns have to exist in the dataframe for plotting
    }
    pred[,thisResp] = 1 #these columns have to exist in the dataframe for plotting
    
    low_color = "#91bfdb"
    high_color = "#fc8d59"
    high_line = "#f46d43"
    low_line = "#4575b4"
    
    if (grad_fill == T) {
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var)) +
        geom_point(alpha = 0.55,stroke = 0, size = point_size) 
      if (!is.null(group_var)) {
        cat("adding lines")
        p1<- p1 + geom_line(aes_string(group = group_var),alpha = .5)
      }
      p1 <- p1 +
        scale_color_gradientn(colors = c(low_line,"grey90",high_line), values = cbar_vals,name = "") +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = "lab"),alpha = .3, linetype = 0) +
        scale_fill_manual(values = c(high_color,low_color)) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",group = "lab"),size = line_size) +
        labs(title = plabels)
    } else {
      black_color = "black"
      green_color = "#4daf4a"
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var)) +
        geom_point(alpha = .35,stroke = 0, size = point_size)
      if (!is.null(group_var)) {
        p1<- p1 + geom_line(aes_string(group = group_var),alpha = .3)
      } 
      p1 <- p1 +
        # scale_color_brewer(type = "qual",palette = "Set1",direction = 1) +
        scale_color_manual(values = c(black_color,green_color)) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = int_var),alpha = .5, linetype = 0) +
        # scale_fill_brewer(type = "qual",palette = "Set1",direction = 1) +
        scale_fill_manual(values = c(black_color,green_color)) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",color = int_var),size = line_size) +
        labs(title = plabels)
    }
  } else {
    
    # No interaction variable, just produce a single line plot
    int_var = "" # This is a lazy solution to making the existing code workable with no int_var.
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    pred <- thisPred %>% select(-init)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    pred[,group_var] = NA
    pred[,thisResp] = 1
    
    p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp)) +
      geom_point(alpha = .3,stroke = 0, size = point_size)
    if (!is.null(group_var)) {
      cat("adding lines")
      p1<- p1 + geom_line(aes_string(group = group_var),alpha = .5)
    }
    p1 <- p1 + geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi"),alpha = .5, linetype = 0) +
      geom_line(data = pred,aes_string(x = smooth_var, y = "fit"),size = line_size) +
      labs(title = plabels)
  }
  
  if (derivative_plot == T) {
    # We will add a bar that shows where the derivative is significant.
    # First make some adjustments to the line plot.
    p1<- p1+theme(text = element_text(size=this_font_size),
                  axis.text = element_text(size = this_font_size),
                  axis.title.y = element_text(size = this_font_size),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.text = element_text(size = this_font_size),
                  legend.title = element_text(size = this_font_size),
                  axis.title = element_text(size = this_font_size),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  plot.margin = unit(c(.2, .2, 0, .2), "cm")) #Top, left,Bottom, right
    scatter <- list(p1)
    
    # Now add the plots using the derivative plotting function
    if (any(grepl(x = row.names(s$s.table),pattern =  ":") & grepl(x=row.names(s$s.table),pattern = int_var))) {
      # Factor levels separately if there is an interaction in the model.
      f<-formula(model) # current formula
      fterms <- terms(f)
      fac <- attr(fterms, "factors")
      idx <- which(as.logical(colSums(fac[grep(x=row.names(fac),pattern = int_var),])))
      new_terms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
      new_formula <- formula(new_terms) # Formula without any interaction terms in the model.
      
      #add derivative gradients for each level of the factor
      num_levels <- length(levels(df[,int_var]))
      level_colors <- suppressWarnings(RColorBrewer::brewer.pal(num_levels,"Set1")) #use the same palette as the line plot
      plotlist = vector(mode = "list",length = num_levels+1) # we will be building a list of plots
      plotlist[1] = scatter # first the scatter plot
      level_colors=c(black_color,green_color)
      for (fcount in 1:num_levels) {
        this_level <- levels(df[,int_var])[fcount]
        df$subset <- df[,int_var] == this_level
        this_mod <- gam(formula = new_formula,data = df,subset = subset)
        # this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        
        if (fcount != num_levels & fcount != 1){
          # get rid of redundant junk
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          this_d$theme$legend.text = element_blank()
          this_d$theme$legend.position="none"
        }
        if (fcount == 1) {
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.text = element_blank()
          this_d$theme$legend.position="none"
        }
        if (fcount == num_levels) {
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          legend <- get_legend(this_d)
        }
        this_d$labels$fill=NULL
        plotlist[fcount+1] = list(this_d+theme(legend.position = "none"))
      }

      pg<-cowplot::plot_grid(rel_heights = c(16*num_levels,rep(num_levels,num_levels-1),2*num_levels),plotlist = plotlist,align = "v",axis = "lr",ncol = 1)
      # final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.1),ncol = 1)
      final_plot <- pg
      print(final_plot)
    } else {
      # No need to split
      d1 <- get_derivs_and_plot(modobj = modobj,smooth_var = smooth_var)
      scatter <- list(p1)
      bar <- list(d1+theme(legend.position = "none"))
      legend <- get_legend(d1+theme(legend.position = "bottom",legend.direction = "horizontal"))
      allplots <- c(scatter,bar)
      pg<-cowplot::plot_grid(rel_heights = c(1,.35),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
      final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.1),ncol = 1)
      print(final_plot)
      final_plot=pg
      
      # print(final_plot)
    }
    
  }    else {
    # No derivative plot
    p1<- p1+theme(text = element_text(size=font_size),
                  axis.text = element_text(size = font_size),
                  legend.text = element_text(size = font_size),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank())
    final_plot <- p1
    print(final_plot)
  }
  
  if (check_diagnostics == T) {
    cp <- check(b,
                a.qq = list(method = "tnorm",
                            a.cipoly = list(fill = "light blue")),
                a.respoi = list(size = 0.5),
                a.hist = list(bins = 10))
    print(cp)
  }
  return(final_plot)
}
