## Use this to start every program.  This clears out previous information from memory
rm(list=ls())

set.seed(121)

## Initalize renv for library lockfile
#library(renv)
#renv::init()

## Packages
#Sys.setenv(RENV_PATHS_RTOOLS = "C:/rtools40/") # https://github.com/rstudio/renv/issues/225

PKG <- c("tidyverse","htmlwidgets","GPareto","patchwork","mco","RColorBrewer")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}
rm(p,PKG)

options(scipen=999) # Prevent scientific notation


# Multi-objective optimization function used for pareto frontiers  --------
# Non-dominated Sorting Genetic Algorithm II using binary based drawing and a repair function to set constraint (penalties and the constraint function don't seem to work)
replifish<-function(n,pops){
  transform_vector_to_exact_ones <- function(v, num_ones) { # Sets constraint by forcing the NSGA2 draws to generate exactly n ones and (D - n) zeros
    # Length of the input vector
    l <- length(v)
    
    # Create a vector of zeros
    transformed_v <- rep(0, l)
    
    # Find indices of the top num_ones values
    top_indices <- order(v, decreasing = TRUE)[1:num_ones]
    
    # Set these positions to 1
    transformed_v[top_indices] <- 1
    
    return(transformed_v)
  }
  
  # z<-runif(nrow(df),0,1) test vector
  
  main.goal1<-function(x){    # x - a vector of indicator variables
    return(sum(transform_vector_to_exact_ones(x,num_ones = n)*df$Total_fish_USD))}
  
  main.goal2<-function(x){    # x - a vector of indicator variables
    return(sum(transform_vector_to_exact_ones(x,num_ones = n)*df$LCOE_MWh*df$`Total MWhyraw`)/sum(transform_vector_to_exact_ones(x,num_ones = n)*df$`Total MWhyraw`))} # Weighted average of LCOE, weighted by energy production
  
  eval<-function(x){
    return(c(main.goal1(x),main.goal2(x)))} # objective function
  
  D<-as.numeric(nrow(df)) # Number of values in the function must equal D
  #n<-63 # Number of cells
  
  system.time(G<-nsga2(fn=eval,
                       idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                       odim=2, # Output dimensions
                       lower.bounds=rep(0,nrow(df)),
                       upper.bounds=rep(1,nrow(df)),
                       popsize=pops,generations=10000, cprob = 0.7, cdist = 5,
                       mprob = 0.2, mdist = 10))
  
  # Identifying pareto sets after optimization
  pareto_indices <- which(G$pareto.optimal == TRUE)
  pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
  names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
  pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
  pareto_solutions$Pareto<-NULL
  pareto_solutions<-unique(pareto_solutions)
  
  
  transform_row_to_exact_ones <- function(row, num_ones) { # Same as vector version above, but for a matrix
    l <- length(row)
    transformed_row <- rep(0, l)
    top_indices <- order(row, decreasing = TRUE)[1:num_ones]
    transformed_row[top_indices] <- 1
    return(transformed_row)
  }
  
  idim_par<-t(apply(G$par, 1, transform_row_to_exact_ones, num_ones = n))
  
  idim_par<-as.data.frame(idim_par[pareto_indices,])
  find_ones_positions <- function(row) {
    return(which(row == 1))
  }
  selected_rows<-t(apply(idim_par,MARGIN = 1, find_ones_positions)) # Lists of selected rows from initial df for pareto sets
  dupes<-as.numeric(nrow(selected_rows)-nrow(unique(selected_rows)))
  selected_rows<-unique(selected_rows) # NSGA2 can produce duplicate pareto sets
  selected_wfgridID<-as.data.frame(selected_rows) %>%
    rowwise() %>%
    mutate(across(everything(), ~df$`Wind farm grid ID`[.x])) %>%
    ungroup()
  # agree<-as.data.frame(table(as.matrix(selected_wfgridID)))
  # agree$Var1<-as.numeric(levels(agree$Var1))[agree$Var1] # https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
  # agree$prob<-agree$Freq/sum(agree$Freq)
  # agree$n<-n
  # agree<-mutate(agree, LCOEMWh = df$LCOE_MWh[match(agree$Var1,df$`Wind farm grid ID`)])
  # agree<-mutate(agree, FishPV = df$Total_fish_USD[match(agree$Var1,df$`Wind farm grid ID`)])
  selected_wfgridID$id<-seq(1,nrow(selected_wfgridID),1) # Adds unique id for each pareto point
  selected_wfgridID<-pivot_longer(data = selected_wfgridID,cols = !id,values_to = "Wind farm grid ID") %>% dplyr::select(!name) # Long format
  
  # Merging fisheries data to pareto points
  dft<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
  dft<-dft %>% dplyr::select(`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD,`Weighted mean LCOE`,`Total MWhyraw`)
  dft<-as.data.frame(dft)
  
  selected_wfgridID<-merge(selected_wfgridID,dft,by="Wind farm grid ID") # Merge
  selected_wfgridID<-unique(selected_wfgridID) # Issue with the duplicate pareto sets command above when n = 1, this catches it
  selected_wfgridID$costofenergy<-selected_wfgridID$`Weighted mean LCOE`*selected_wfgridID$`Total MWhyraw`
  LCOEMWh<-sum(selected_wfgridID$costofenergy)/sum(selected_wfgridID$`Total MWhyraw`)*1000 # Mean LCOE for the development scenario (energy production weighted)
  s_ids<-selected_wfgridID # Saves selected pareto sets with ids to a dataframe for return
  s_ids$n<-n
  
  selected_wfgridID$`Weighted mean LCOE`<-NULL
  selected_wfgridID$`Total MWhyraw`<-NULL
  selected_wfgridID$costofenergy<-NULL
  selected_wfgridID<-selected_wfgridID %>% dplyr::select(!"Wind farm grid ID") %>% group_by(id) %>% # Sums exposure per species for a given pareto set 
    summarise(across(everything(), list(sum),.names = "sum_{.col}"))
  selected_wfgridIDmn<-selected_wfgridID %>% dplyr::select(!id) %>% # Mean exposure per species for all pareto sets
    summarise(across(everything(), list(mean),.names = "mean_{.col}")) 
  selected_wfgridIDmd<-selected_wfgridID %>% dplyr::select(!id) %>% # Median exposure per species for all pareto sets
    summarise(across(everything(), list(median),.names = "median_{.col}")) 
  
  # Function to calculate confidence interval using t.test
  calculate_ci <- function(column) {
    t_result <- t.test(column)
    lower <- t_result$conf.int[1]
    upper <- t_result$conf.int[2]
    
    # Check for NaN and replace with mean for when there is only a single value in the entire sample
    if (is.nan(lower)) lower <- as.numeric(t_result$estimate)
    if (is.nan(upper)) upper <- as.numeric(t_result$estimate)
    
    return(c(lower, upper))
  }
  ci_df<-as.data.frame(t(sapply(selected_wfgridID[, setdiff(names(selected_wfgridID), "id")], calculate_ci)))
  colnames(ci_df)<-c("Lower", "Upper")
  ci_df$Fishery<-rownames(ci_df)
  ci_df$Fishery<-gsub("sum_", "",ci_df$Fishery)
  
  bigfishmn<-selected_wfgridIDmn %>% pivot_longer(cols = everything(),names_to = "Fishery",values_to = "Mean PV", names_prefix = "mean_sum_")
  bigfishmd<-selected_wfgridIDmd %>% pivot_longer(cols = everything(),names_to = "Fishery",values_to = "Median PV", names_prefix = "median_sum_")
  bigfish<-merge(bigfishmn,bigfishmd, by = "Fishery")
  bigfish<-merge(bigfish,ci_df,by="Fishery")
  
  # m<-15000 # Highest value for expected values of n
  # bigfish<-data.frame("Dungeness_USD"=rep(NA,m),"At-sea_hake_USD"=rep(NA,m),"Shore_hake_USD"=rep(NA,m),"Market_squid_USD"=rep(NA,m),"Pink_shrimp_USD"=rep(NA,m),"Albacore_USD"=rep(NA,m),"Chinook_USD"=rep(NA,m),"Sablefish_USD"=rep(NA,m),"Spiny_lobster_USD"=rep(NA,m),"LCOEMWh"=rep(NA,m),"Sites"=rep(NA,m))
  # for(j in 1:m){ # Sums species PV for a drawn sample of size n and iterates across samples while tracking draws in "Sites" column
  #   s_agree<-sample(agree$Var1,n,replace = FALSE,prob = agree$prob)
  #   fishroll<-data.frame("Dungeness_USD"=rep(NA,n),"At-sea_hake_USD"=rep(NA,n),"Shore_hake_USD"=rep(NA,n),"Market_squid_USD"=rep(NA,n),"Pink_shrimp_USD"=rep(NA,n),"Albacore_USD"=rep(NA,n),"Chinook_USD"=rep(NA,n),"Sablefish_USD"=rep(NA,n),"Spiny_lobster_USD"=rep(NA,n),"Weighted mean LCOE"=rep(NA,n))
  #   for(i in 1:n){ # Identifies species PV from a given weighted sample of size n drawn from agreement vector
  #     fishroll[i,]<-dft[dft$`Wind farm grid ID`==s_agree[i],2:11]
  #   }
  #   LCOEMWh<-as.data.frame(mean(fishroll$Weighted.mean.LCOE)*1000)
  #   fishroll$Weighted.mean.LCOE<-NULL
  #   fishroll<-as.data.frame(t(colSums(fishroll))) #paste(s_agree, collapse = " ")
  #   fishroll<-cbind(fishroll,LCOEMWh,paste(s_agree, collapse = " "))
  #   names(fishroll)[names(fishroll) == 'paste(s_agree, collapse = " ")']<-"Sites"
  #   names(fishroll)[names(fishroll) == 'mean(fishroll$Weighted.mean.LCOE) * 1000']<-"LCOEMWh"
  #   bigfish[j,]<-fishroll
  # }
  # 
  # LCOEMWh<-mean(bigfish$LCOEMWh)
  # bigfish$LCOEMWh<-NULL
  # bigfish<-bigfish %>% pivot_longer(cols = !Sites,names_to = "Fishery",values_to = "PV")
  # bigfish$Fishery<-gsub("_USD$", "",bigfish$Fishery)
  # bigfish$Fishery<-gsub("_", " ",bigfish$Fishery)
  # bigfish$Fishery<-gsub("\\.", "-",bigfish$Fishery)
  # 
  # bigfish<-bigfish %>% group_by(Fishery) %>% 
  #   summarise(meanfishPV = mean(PV))
  
  bigfish$n<-n
  bigfish$`Mean PV`<-ifelse(bigfish$n==1,bigfish$`Mean PV`/9,bigfish$`Mean PV`) # Handles formatting issue with n == 1
  bigfish$`Median PV`<-ifelse(bigfish$n==1,bigfish$`Median PV`/9,bigfish$`Median PV`) # Handles formatting issue with n == 1
  bigfish$LCOEMWh<-LCOEMWh
  bigfish$dupes<-dupes
  
  return(list(bigfish = bigfish,s_ids = s_ids,pareto_solutions = pareto_solutions))
  #return(list(bigfish = bigfish))
}

# Wind results (supports numbers in results text) ------------------------------------------------------------
df<-read_csv("OWEP LCOE NPV & fishing PV data V5 DO NOT DISTRIBUTE.csv")
df<-df[df$LCOE_kWh>0,] %>% drop_na(LCOE_kWh)
df$LCOE_MWh<-df$LCOE_kWh*1000
df$LCOE_kWh<-NULL
mean(df$LCOE_MWh)
df %>% group_by(State) %>% 
  summarise(mean = mean(LCOE_MWh), std = sd(LCOE_MWh))


# Scatterplot and table of LCOE per fishery Fig 2 -------------------------
df<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
df<-df %>% drop_na(`Mean LCOE_kWh`)
df<-df[df$`OWEP grid cells (n)`>34,]
df$Total_fish_USD<-NULL
df<-pivot_longer(df,cols = ends_with("_USD"),names_to = "fishery",values_to = "PV")
df<-df[df$PV>0,]
#df$Total_fish_USD<-df$Total_fish_USD/1000000
df$LCOE_MWh<-df$`Weighted mean LCOE`*1000
df<-df %>% dplyr::select(`Wind farm grid ID`,`Wind farm state`,LCOE_MWh,fishery,PV)
df<-as.data.frame(df)
df$`Wind farm state`<-as.factor(df$`Wind farm state`)
df$fishery<-gsub("_USD$", "",df$fishery)
df$fishery<-gsub("_", " ",df$fishery)
#df$fishery<-as.factor(df$fishery)

sorted_fish<-df %>% distinct(fishery) %>% 
  arrange(desc(fishery)) %>% unlist()

df$fishery<-factor(df$fishery, levels = sorted_fish)

ggplot() + # Figure 2
  #geom_point(data = df,aes(x = PV, y = LCOE_MWh, color = fishery)) + # Scatterplot of fisheries PV and LCOE is too jumbled with extreme values
  geom_boxplot(data = df, aes(y = fishery, x = LCOE_MWh),outlier.shape = NA) +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)", color = "Fishery") +
  labs(y="",x = "LCOE ($/MWh)") + xlim(c(50,500)) +
  theme_minimal()


# Pareto frontiers for 2030 and 2045 targets ------------------------------
# Loading data
df<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
df<-df %>%
  dplyr::select(`Wind farm grid ID`,`Wind farm state`,`OWEP grid cells (n)`,Total_fish_USD,`Weighted mean LCOE`,`Total MWhyraw`,`Wind farm area (m2)`)
df<-df %>% drop_na(`Weighted mean LCOE`)
df<-df[df$`OWEP grid cells (n)`>34,] # Choosing 34 as a cutoff for what could be considered full cells, could amend
df$Total_fish_USD<-df$Total_fish_USD/1000000
df$LCOE_MWh<-df$`Weighted mean LCOE`*1000
df$`Weighted mean LCOE`<-NULL
df$`Wind farm area (m2)`<-NULL # This isn't a correct calculation, this is actually cell size
df<-as.data.frame(df)
df2<-df

# Targets: 
# CA is 5GW by 2030 and 25GW by 2045, OR is 3GW by 2030. Could simulate WA at 3GW. and 
# For 2030, assuming 3MW per km2 and each cell is 295km2 (36*8.185km2 cells): CA is 1667km2 (6 cells), OR is 1000km2 (4 cells), WA is 1000km2 (4 cells), and together is 3667km2 (13 cells)
# For 2045, assuming a 5 times increase in GW targets: Region 63 cells, CA 28 cells, OR 17 cells, WA 17 cells

# Region 2030
df<-df2
ns<-13
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))
 
pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))
    
s_idsR<-s_ids  
pareto_solutionsR<-pareto_solutions
biggerfishR<-do.call(rbind, lapply(biggerfish, `[[`, 1))

region<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none") +
  coord_cartesian(xlim = c(0, 80))

# CA 2030
df<-df2[df2$`Wind farm state`=="CA",]
ns<-6
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))

pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))

s_idsCA<-s_ids  
pareto_solutionsCA<-pareto_solutions
biggerfishCA<-do.call(rbind, lapply(biggerfish, `[[`, 1))

CA<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  theme_minimal() +
  theme(legend.title = element_blank(),legend.position="none") +
  coord_cartesian(xlim = c(0, 8))

# OR 2030
df<-df2[df2$`Wind farm state`=="OR",]
ns<-4
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))

pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))

s_idsOR<-s_ids  
pareto_solutionsOR<-pareto_solutions
biggerfishOR<-do.call(rbind, lapply(biggerfish, `[[`, 1))

OR<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  theme_minimal() +
  theme(legend.title = element_blank(),legend.position="none")

# WA 2030
df<-df2[df2$`Wind farm state`=="WA",]
ns<-4
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))

pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))

s_idsWA<-s_ids  
pareto_solutionsWA<-pareto_solutions
biggerfishWA<-do.call(rbind, lapply(biggerfish, `[[`, 1))

WA<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  theme_minimal() +
  theme(legend.title = element_blank(),legend.position="none")

# Region 2045
df<-df2
ns<-63
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))

pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))

s_idsR45<-s_ids  
pareto_solutionsR45<-pareto_solutions
biggerfishR45<-do.call(rbind, lapply(biggerfish, `[[`, 1))

region45<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  theme_minimal() +
  theme(legend.title = element_blank(),legend.position="none")  
  
# CA 2045
df<-df2[df2$`Wind farm state`=="CA",]
ns<-28
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))

pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))

s_idsCA45<-s_ids  
pareto_solutionsCA45<-pareto_solutions
biggerfishCA45<-do.call(rbind, lapply(biggerfish, `[[`, 1))

CA45<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  theme_minimal() +
  theme(legend.title = element_blank(),legend.position="none")

# OR 2045
df<-df2[df2$`Wind farm state`=="OR",]
ns<-17
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))

pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))

s_idsOR45<-s_ids  
pareto_solutionsOR45<-pareto_solutions
biggerfishOR45<-do.call(rbind, lapply(biggerfish, `[[`, 1))

OR45<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  theme_minimal() +
  theme(legend.title = element_blank(),legend.position="none")

# WA 2045
df<-df2[df2$`Wind farm state`=="WA",]
ns<-17
pop<-1000
system.time(biggerfish<-map2(ns,pop,replifish))

pareto_solutions<-do.call(rbind, lapply(biggerfish, `[[`, 3))
pareto_solutions_ordered<-pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ] # Prep for plotting
s_ids<-do.call(rbind, lapply(biggerfish, `[[`, 2))

s_ids<-s_ids %>% 
  mutate(TotRev = rowSums(select(., ends_with("_USD")))) %>% # Adding up columns that end in "_USD" only - flexible to some columns not being present
  group_by(id) %>% 
  summarise(`Sum fishing PV ($Mil)` = sum(TotRev)/1000000, `Mean LCOE ($/MWh)` = sum(costofenergy)/sum(`Total MWhyraw`)*1000, GridId = paste(`Wind farm grid ID`, collapse = ", ")) %>% 
  mutate(rank = rank(`Sum fishing PV ($Mil)`)) %>% 
  filter(rank %in% c(round(length(rank)*.1),round(length(rank)*.5),round(length(rank)*.9))) %>%   # Keeping the 10th, 50th, and 90th percentile observations
  mutate(lab = ifelse(`Sum fishing PV ($Mil)`==max(`Sum fishing PV ($Mil)`),"Wind",
                      ifelse(`Sum fishing PV ($Mil)`==min(`Sum fishing PV ($Mil)`),"Fish","Balanced")))

s_idsWA45<-s_ids  
pareto_solutionsWA45<-pareto_solutions
biggerfishWA45<-do.call(rbind, lapply(biggerfish, `[[`, 1))

WA45<-ggplot() +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  geom_point(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), size = 5) +
  geom_text(data = s_ids, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, label = lab), hjust = -0.2, vjust = -.75) +
  theme_minimal() +
  theme(legend.title = element_blank(),legend.position="none")

# Grouping PPF plots
#region + CA + OR + WA

region<-region + labs(x="Sum fishing present value exposure ($ Mil)", title = "Region 2030") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot") # Relabeling x axis
region45<-region45 + labs(x="Sum fishing present value exposure ($ Mil)", title = "Region 2045") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot")
CA<-CA + labs(x="Sum fishing present value exposure ($ Mil)", title = "CA 2030") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot") + coord_cartesian(xlim = c(0, 8))
CA45<-CA45 + labs(x="Sum fishing present value exposure ($ Mil)", title = "CA 2045") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot")
WA<-WA + labs(x="Sum fishing present value exposure ($ Mil)", title = "OR 2030") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot") + coord_cartesian(xlim = c(0, 350))
WA45<-WA45 + labs(x="Sum fishing present value exposure ($ Mil)", title = "OR 2045") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot")
OR<-OR + labs(x="Sum fishing present value exposure ($ Mil)", title = "WA 2030") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot")
OR45<-OR45 + labs(x="Sum fishing present value exposure ($ Mil)", title = "WA 2045") + theme(plot.title = element_text(hjust = 0.5,vjust = -1.5), plot.title.position = "plot")

#region + region + region45 + CA45 + OR45 + WA45 + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect")
#region + region45 + CA45 + OR45 + WA45 + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect")
#region + region45 + CA + OR + WA + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect", design = "AAABBB\nCCDDEE")

# a1<-region + region45 + plot_layout(axis_titles = "collect")
# a2<-CA + OR + WA + plot_layout(axis_titles = "collect")
# a1 / a2 + plot_annotation(tag_levels = list(c("A","B","C","D","E"))) # Old fig 4
# 
# region45 + CA45 + OR45 + WA45 + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect") # Old supplement figure

region + region45 + CA + CA45 + OR + OR45 + WA + WA45 + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect", ncol = 2) # New figure 4 that combines old main and supplement figures and adds points for mapping

#region + region + region45 + CA + OR + WA + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect")
#(region + region45) / (CA + OR + WA) + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect")

bigfish30 / bigfish45 + plot_annotation(tag_levels = 'A')

bigfish30perc / bigfish45perc + plot_annotation(tag_levels = 'A')

# Table of wind cell ids for mapping results represnted by the points on the frontiers ala https://doi.org/10.1016/j.biocon.2008.03.022 Fig. 3
s_idsCA$area<-"CA"
s_idsCA45$area<-"CA"
s_idsOR$area<-"OR"
s_idsOR45$area<-"OR"
s_idsWA$area<-"WA"
s_idsWA45$area<-"WA"
s_idsR$area<-"Region"
s_idsR45$area<-"Region"

s_idsCA$target<-2030
s_idsCA45$target<-2045
s_idsOR$target<-2030
s_idsOR45$target<-2045
s_idsWA$target<-2030
s_idsWA45$target<-2045
s_idsR$target<-2030
s_idsR45$target<-2045

s_ids<-rbind(s_idsCA,s_idsCA45,s_idsOR,s_idsOR45,s_idsR,s_idsR45,s_idsWA,s_idsWA45)
s_ids$id<-NULL

s_ids<-s_ids %>% separate_rows(GridId, sep = ",") # Long format
s_ids$GridId<-as.numeric(s_ids$GridId)

write.csv(s_ids,"Pareto_mapping_ids.csv")


# Fisheries exposure replicated across many targets -----------------------
## Region
# Loading data
df<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
df<-df %>%
  dplyr::select(`Wind farm grid ID`,`Wind farm state`,`OWEP grid cells (n)`,Total_fish_USD,`Weighted mean LCOE`,`Total MWhyraw`,`Wind farm area (m2)`)
df<-df %>% drop_na(`Weighted mean LCOE`)
df<-df[df$`OWEP grid cells (n)`>34,] # Choosing 34 as a cutoff for what could be considered full cells, could amend
df$Total_fish_USD<-df$Total_fish_USD/1000000
df$LCOE_MWh<-df$`Weighted mean LCOE`*1000
df$`Weighted mean LCOE`<-NULL
df$`Wind farm area (m2)`<-NULL # This isn't a correct calculation, this is actually cell size
df<-as.data.frame(df)

ns<-c(seq(2,62,10))
pop<-rep(1000,length(ns))

#system.time(biggerfish<-map2_dfr(ns,pop,replifish)) # Use if the replifish function returns a single dataframe
system.time(biggerfish<-map2(ns,pop,replifish)) # Use if the replifish function returns a list

s_idsR<-do.call(rbind, lapply(biggerfish, `[[`, 2))
biggerfish<-do.call(rbind, lapply(biggerfish, `[[`, 1))

# biggerfish %>% group_by(n) %>% 
#   summarise(SumMean = sum(`Mean PV`), SumMedian = sum(`Median PV`), MeanLCOE = mean(LCOEMWh)) %>% 
#   print(n = 100)

palette <- brewer.pal(n = length(unique(biggerfish$Fishery)), name = "Paired")

dft<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
dft<-dft %>% dplyr::select(`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD,`Weighted mean LCOE`)
dft<-as.data.frame(dft)
fishsum<-dft %>% pivot_longer(cols = !`Wind farm grid ID`,names_to = "Fishery",values_to = "PV") %>% 
  group_by(Fishery) %>% 
  summarise(fishsum = sum(PV))

fishsum$Fishery<-gsub("_USD$", "",fishsum$Fishery)
fishsum$Fishery<-gsub("_", " ",fishsum$Fishery)
fishsum$Fishery<-gsub("\\.", "-",fishsum$Fishery)

biggerfish$Fishery<-gsub("_USD$", "",biggerfish$Fishery)
biggerfish$Fishery<-gsub("_", " ",biggerfish$Fishery)
biggerfish$Fishery<-gsub("\\.", "-",biggerfish$Fishery)

biggerfish<-merge(biggerfish,fishsum,by="Fishery")
biggerfish$PVmnperc<-(biggerfish$`Mean PV`/biggerfish$fishsum)*100
biggerfish$PVmdperc<-(biggerfish$`Median PV`/biggerfish$fishsum)*100
biggerfish$PVlperc<-(biggerfish$Lower/biggerfish$fishsum)*100
biggerfish$PVuperc<-(biggerfish$Upper/biggerfish$fishsum)*100

biggerfishR<-biggerfish

bf_nom_mn<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,62))

bf_nom_md<-ggplot(data = biggerfish, aes(x=n*.9, y=`Median PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median present value ($Mil)", colour = "Fishery") +
  xlim(c(0,62))

l_nom<-ggplot(data = biggerfish, aes(x=n*.9, y=LCOEMWh)) + # LCOE expected impact across targets, weird artifacts from lots of pareto sites with high LCOE (TRY MEDIAN)
  geom_line(linewidth = 1)

bf_perc_mn<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,62))

bf_perc_md<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmdperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median Percent", colour = "Fishery") +
  xlim(c(0,62))

ci_bf_nom_mn<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Lower/1000000, ymax = Upper/1000000),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,62))

ci_bf_perc_mn<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = PVlperc, ymax = PVuperc),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,62))

## CA
# Loading data
df<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
df<-df %>%
  dplyr::select(`Wind farm grid ID`,`Wind farm state`,`OWEP grid cells (n)`,Total_fish_USD,`Weighted mean LCOE`,`Total MWhyraw`,`Wind farm area (m2)`)
df<-df %>% drop_na(`Weighted mean LCOE`)
df<-df[df$`OWEP grid cells (n)`>34,] # Choosing 34 as a cutoff for what could be considered full cells, could amend
df$Total_fish_USD<-df$Total_fish_USD/1000000
df$LCOE_MWh<-df$`Weighted mean LCOE`*1000
df$`Weighted mean LCOE`<-NULL
df$`Wind farm area (m2)`<-NULL # This isn't a correct calculation, this is actually cell size
df<-as.data.frame(df)
df<-df %>% filter(`Wind farm state`=="CA")

ns<-c(seq(5,30,5))
pop<-rep(1000,length(ns))

#system.time(biggerfish<-map2_dfr(ns,pop,replifish)) # Use if the replifish function returns a single dataframe
system.time(biggerfish<-map2(ns,pop,replifish)) # Use if the replifish function returns a list

s_idsCA<-do.call(rbind, lapply(biggerfish, `[[`, 2))
biggerfish<-do.call(rbind, lapply(biggerfish, `[[`, 1))

biggerfish %>% group_by(n) %>% 
  summarise(SumMean = sum(`Mean PV`), SumMedian = sum(`Median PV`), MeanLCOE = mean(LCOEMWh)) %>%
  print(n = 100)

palette <- brewer.pal(n = length(unique(biggerfish$Fishery)), name = "Paired")

dft<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
dft<-dft %>% dplyr::select(`Wind farm state`,`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD,`Weighted mean LCOE`)
dft<-as.data.frame(dft)
dft<-dft %>% filter(`Wind farm state`=="CA") %>% select(!`Wind farm state`)
fishsum<-dft %>% pivot_longer(cols = !`Wind farm grid ID`,names_to = "Fishery",values_to = "PV") %>% 
  group_by(Fishery) %>% 
  summarise(fishsum = sum(PV))

fishsum$Fishery<-gsub("_USD$", "",fishsum$Fishery)
fishsum$Fishery<-gsub("_", " ",fishsum$Fishery)
fishsum$Fishery<-gsub("\\.", "-",fishsum$Fishery)

biggerfish$Fishery<-gsub("_USD$", "",biggerfish$Fishery)
biggerfish$Fishery<-gsub("_", " ",biggerfish$Fishery)
biggerfish$Fishery<-gsub("\\.", "-",biggerfish$Fishery)

biggerfish<-merge(biggerfish,fishsum,by="Fishery")
biggerfish$PVmnperc<-(biggerfish$`Mean PV`/biggerfish$fishsum)*100
biggerfish$PVmdperc<-(biggerfish$`Median PV`/biggerfish$fishsum)*100
biggerfish$PVlperc<-(biggerfish$Lower/biggerfish$fishsum)*100
biggerfish$PVuperc<-(biggerfish$Upper/biggerfish$fishsum)*100
biggerfish<-biggerfish %>% mutate(across(everything(), ~replace(.x, is.nan(.x), 0))) # Some species with no harvest in state waters introduces a NaN when dividing 0 by 0

biggerfishCA<-biggerfish

bf_nom_mnCA<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 5, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 25, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 25, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,30))

bf_nom_mdCA<-ggplot(data = biggerfish, aes(x=n*.9, y=`Median PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 5, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 25, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 25, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median present value ($Mil)", colour = "Fishery") +
  xlim(c(0,30))

l_nomCA<-ggplot(data = biggerfish, aes(x=n*.9, y=LCOEMWh)) + # LCOE expected impact across targets, weird artifacts from lots of pareto sites with high LCOE (TRY MEDIAN)
  geom_line(linewidth = 1)

bf_perc_mnCA<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 5, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 25, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 25, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,30))

bf_perc_mdCA<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmdperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 5, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 25, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 25, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median Percent", colour = "Fishery") +
  xlim(c(0,30))

ci_bf_nom_mnCA<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Lower/1000000, ymax = Upper/1000000),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 5, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 25, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 25, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,30))

ci_bf_perc_mnCA<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = PVlperc, ymax = PVuperc),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 5, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 25, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 25, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,30))

## OR
# Loading data
df<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
df<-df %>%
  dplyr::select(`Wind farm grid ID`,`Wind farm state`,`OWEP grid cells (n)`,Total_fish_USD,`Weighted mean LCOE`,`Total MWhyraw`,`Wind farm area (m2)`)
df<-df %>% drop_na(`Weighted mean LCOE`)
df<-df[df$`OWEP grid cells (n)`>34,] # Choosing 34 as a cutoff for what could be considered full cells, could amend
df$Total_fish_USD<-df$Total_fish_USD/1000000
df$LCOE_MWh<-df$`Weighted mean LCOE`*1000
df$`Weighted mean LCOE`<-NULL
df$`Wind farm area (m2)`<-NULL # This isn't a correct calculation, this is actually cell size
df<-as.data.frame(df)
df<-df %>% filter(`Wind farm state`=="OR")

ns<-c(seq(3,18,3))
pop<-rep(1000,length(ns))

#system.time(biggerfish<-map2_dfr(ns,pop,replifish)) # Use if the replifish function returns a single dataframe
system.time(biggerfish<-map2(ns,pop,replifish)) # Use if the replifish function returns a list

s_idsOR<-do.call(rbind, lapply(biggerfish, `[[`, 2))
biggerfish<-do.call(rbind, lapply(biggerfish, `[[`, 1))

biggerfish %>% group_by(n) %>% 
  summarise(SumMean = sum(`Mean PV`), SumMedian = sum(`Median PV`), MeanLCOE = mean(LCOEMWh)) %>% 
  print(n = 100)

palette <- brewer.pal(n = length(unique(biggerfish$Fishery)), name = "Paired")

dft<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
dft<-dft %>% dplyr::select(`Wind farm state`,`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD,`Weighted mean LCOE`)
dft<-as.data.frame(dft)
dft<-dft %>% filter(`Wind farm state`=="OR") %>% select(!`Wind farm state`)
fishsum<-dft %>% pivot_longer(cols = !`Wind farm grid ID`,names_to = "Fishery",values_to = "PV") %>% 
  group_by(Fishery) %>% 
  summarise(fishsum = sum(PV))

fishsum$Fishery<-gsub("_USD$", "",fishsum$Fishery)
fishsum$Fishery<-gsub("_", " ",fishsum$Fishery)
fishsum$Fishery<-gsub("\\.", "-",fishsum$Fishery)

biggerfish$Fishery<-gsub("_USD$", "",biggerfish$Fishery)
biggerfish$Fishery<-gsub("_", " ",biggerfish$Fishery)
biggerfish$Fishery<-gsub("\\.", "-",biggerfish$Fishery)

biggerfish<-merge(biggerfish,fishsum,by="Fishery")
biggerfish$PVmnperc<-(biggerfish$`Mean PV`/biggerfish$fishsum)*100
biggerfish$PVmdperc<-(biggerfish$`Median PV`/biggerfish$fishsum)*100
biggerfish$PVlperc<-(biggerfish$Lower/biggerfish$fishsum)*100
biggerfish$PVuperc<-(biggerfish$Upper/biggerfish$fishsum)*100
biggerfish<-biggerfish %>% mutate(across(everything(), ~replace(.x, is.nan(.x), 0))) # Some species with no harvest in state waters introduces a NaN when dividing 0 by 0

biggerfishOR<-biggerfish

bf_nom_mnOR<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,18))

bf_nom_mdOR<-ggplot(data = biggerfish, aes(x=n*.9, y=`Median PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median present value ($Mil)", colour = "Fishery") +
  xlim(c(0,18))

l_nomOR<-ggplot(data = biggerfish, aes(x=n*.9, y=LCOEMWh)) + # LCOE expected impact across targets, weird artifacts from lots of pareto sites with high LCOE (TRY MEDIAN)
  geom_line(linewidth = 1)

bf_perc_mnOR<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,18))

bf_perc_mdOR<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmdperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median Percent", colour = "Fishery") +
  xlim(c(0,18))

ci_bf_nom_mnOR<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Lower/1000000, ymax = Upper/1000000),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,18))

ci_bf_perc_mnOR<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = PVlperc, ymax = PVuperc),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,18))

## WA
# Loading data
df<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
df<-df %>%
  dplyr::select(`Wind farm grid ID`,`Wind farm state`,`OWEP grid cells (n)`,Total_fish_USD,`Weighted mean LCOE`,`Total MWhyraw`,`Wind farm area (m2)`)
df<-df %>% drop_na(`Weighted mean LCOE`)
df<-df[df$`OWEP grid cells (n)`>34,] # Choosing 34 as a cutoff for what could be considered full cells, could amend
df$Total_fish_USD<-df$Total_fish_USD/1000000
df$LCOE_MWh<-df$`Weighted mean LCOE`*1000
df$`Weighted mean LCOE`<-NULL
df$`Wind farm area (m2)`<-NULL # This isn't a correct calculation, this is actually cell size
df<-as.data.frame(df)
df<-df %>% filter(`Wind farm state`=="WA")

ns<-c(seq(3,18,3))
pop<-rep(1000,length(ns))

#system.time(biggerfish<-map2_dfr(ns,pop,replifish)) # Use if the replifish function returns a single dataframe
system.time(biggerfish<-map2(ns,pop,replifish)) # Use if the replifish function returns a list

s_idsWA<-do.call(rbind, lapply(biggerfish, `[[`, 2))
biggerfish<-do.call(rbind, lapply(biggerfish, `[[`, 1))

biggerfish %>% group_by(n) %>% 
  summarise(SumMean = sum(`Mean PV`), SumMedian = sum(`Median PV`), MeanLCOE = mean(LCOEMWh)) %>% 
  print(n = 100)

palette <- brewer.pal(n = length(unique(biggerfish$Fishery)), name = "Paired")

dft<-read_csv("OWEP output & fishing PV data V5 DO NOT DISTRIBUTE.csv")
dft<-dft %>% dplyr::select(`Wind farm state`,`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD,`Weighted mean LCOE`)
dft<-as.data.frame(dft)
dft<-dft %>% filter(`Wind farm state`=="WA") %>% select(!`Wind farm state`)
fishsum<-dft %>% pivot_longer(cols = !`Wind farm grid ID`,names_to = "Fishery",values_to = "PV") %>% 
  group_by(Fishery) %>% 
  summarise(fishsum = sum(PV))

fishsum$Fishery<-gsub("_USD$", "",fishsum$Fishery)
fishsum$Fishery<-gsub("_", " ",fishsum$Fishery)
fishsum$Fishery<-gsub("\\.", "-",fishsum$Fishery)

biggerfish$Fishery<-gsub("_USD$", "",biggerfish$Fishery)
biggerfish$Fishery<-gsub("_", " ",biggerfish$Fishery)
biggerfish$Fishery<-gsub("\\.", "-",biggerfish$Fishery)

biggerfish<-merge(biggerfish,fishsum,by="Fishery")
biggerfish$PVmnperc<-(biggerfish$`Mean PV`/biggerfish$fishsum)*100
biggerfish$PVmdperc<-(biggerfish$`Median PV`/biggerfish$fishsum)*100
biggerfish$PVlperc<-(biggerfish$Lower/biggerfish$fishsum)*100
biggerfish$PVuperc<-(biggerfish$Upper/biggerfish$fishsum)*100
biggerfish<-biggerfish %>% mutate(across(everything(), ~replace(.x, is.nan(.x), 0))) # Some species with no harvest in state waters introduces a NaN when dividing 0 by 0

biggerfishWA<-biggerfish

bf_nom_mnWA<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,18))

bf_nom_mdWA<-ggplot(data = biggerfish, aes(x=n*.9, y=`Median PV`/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(`Median PV`)/1000000 - diff(range(`Median PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median present value ($Mil)", colour = "Fishery") +
  xlim(c(0,18))

l_nomWA<-ggplot(data = biggerfish, aes(x=n*.9, y=LCOEMWh)) + # LCOE expected impact across targets, weird artifacts from lots of pareto sites with high LCOE (TRY MEDIAN)
  geom_line(linewidth = 1)

bf_perc_mnWA<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,18))

bf_perc_mdWA<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmdperc, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(PVmdperc) - diff(range(PVmdperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Median Percent", colour = "Fishery") +
  xlim(c(0,18))

ci_bf_nom_mnWA<-ggplot(data = biggerfish, aes(x=n*.9, y=`Mean PV`/1000000, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Lower/1000000, ymax = Upper/1000000),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(`Mean PV`)/1000000 - diff(range(`Mean PV`)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,18))

ci_bf_perc_mnWA<-ggplot(data = biggerfish, aes(x=n*.9, y=PVmnperc, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = PVlperc, ymax = PVuperc),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 3, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 15, y = min(PVmnperc) - diff(range(PVmnperc)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,18))

## Gathering figures
# Combined figure - adds up exposure when forced to develop to targets within state-adjacent waters
biggerfishCA$state<-"CA"
biggerfishOR$state<-"OR"
biggerfishWA$state<-"WA"
biggerfishcomb<-rbind(biggerfishCA,biggerfishOR,biggerfishWA)
# Grouping factor
biggerfishcomb$id<-ifelse(biggerfishcomb$n==3&biggerfishcomb$state %in% c("OR","WA") | biggerfishcomb$n==5&biggerfishcomb$state %in% c("CA"),1,
                          ifelse(biggerfishcomb$n==6&biggerfishcomb$state %in% c("OR","WA") | biggerfishcomb$n==10&biggerfishcomb$state %in% c("CA"),2,
                                 ifelse(biggerfishcomb$n==9&biggerfishcomb$state %in% c("OR","WA") | biggerfishcomb$n==15&biggerfishcomb$state %in% c("CA"),3,
                                        ifelse(biggerfishcomb$n==12&biggerfishcomb$state %in% c("OR","WA") | biggerfishcomb$n==20&biggerfishcomb$state %in% c("CA"),4,
                                               ifelse(biggerfishcomb$n==15&biggerfishcomb$state %in% c("OR","WA") | biggerfishcomb$n==25&biggerfishcomb$state %in% c("CA"),5,
                                                      ifelse(biggerfishcomb$n==18&biggerfishcomb$state %in% c("OR","WA") | biggerfishcomb$n==30&biggerfishcomb$state %in% c("CA"),6,NA))))))

biggerfishcomb$CIradius<-biggerfishcomb$`Mean PV` - biggerfishcomb$Lower # https://stats.stackexchange.com/questions/223924/how-to-add-up-partial-confidence-intervals-to-create-a-total-confidence-interval

# Sum by group
bf_comb<-biggerfishcomb %>% 
  group_by(id,Fishery) %>% 
  summarise(PVnomsum_mn = sum(`Mean PV`),PVnomsum_md = (sum(`Median PV`)^.5), PVCI = (sum(CIradius^2))^.5, GW = sum(n)*.9) %>% 
  as.data.frame()

# Have percentages be expressed as percent of regional fishery value, not summed state percentages
rfv<-biggerfishR %>% dplyr::select(Fishery,fishsum) %>% distinct()
bf_comb<-merge(bf_comb,rfv,by = "Fishery")
bf_comb$PVpercsum_mn<-(bf_comb$PVnomsum_mn/bf_comb$fishsum)*100
bf_comb$PVpercsum_md<-(bf_comb$PVnomsum_md/bf_comb$fishsum)*100

# Confidence intervals
bf_comb$Lower<-bf_comb$PVnomsum_mn-bf_comb$PVCI
bf_comb$Upper<-bf_comb$PVnomsum_mn+bf_comb$PVCI
bf_comb$PVlperc<-bf_comb$PVpercsum_mn-((bf_comb$PVCI/bf_comb$fishsum)*100)
bf_comb$PVuperc<-bf_comb$PVpercsum_mn+((bf_comb$PVCI/bf_comb$fishsum)*100)

# Plots
bf_nom_mns<-ggplot(data = bf_comb, aes(x=GW, y=PVnomsum_mn/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVnomsum_mn)/1000000 - diff(range(PVnomsum_mn)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVnomsum_mn)/1000000 - diff(range(PVnomsum_mn)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,60))

bf_nom_mds<-ggplot(data = bf_comb, aes(x=GW, y=PVnomsum_md/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVnomsum_md)/1000000 - diff(range(PVnomsum_md)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVnomsum_md)/1000000 - diff(range(PVnomsum_md)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,60))

bf_perc_mns<-ggplot(data = bf_comb, aes(x=GW, y=PVpercsum_mn, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVpercsum_mn) - diff(range(PVpercsum_mn)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVpercsum_mn) - diff(range(PVpercsum_mn)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,60))

bf_perc_mds<-ggplot(data = bf_comb, aes(x=GW, y=PVpercsum_md, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), linewidth = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVpercsum_md) - diff(range(PVpercsum_md)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVpercsum_md) - diff(range(PVpercsum_md)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,60))

ci_bf_nom_mns<-ggplot(data = bf_comb, aes(x=GW, y=PVnomsum_mn/1000000, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Lower/1000000, ymax = Upper/1000000),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVnomsum_mn)/1000000 - diff(range(PVnomsum_mn)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVnomsum_mn)/1000000 - diff(range(PVnomsum_mn)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery") +
  xlim(c(0,60))

ci_bf_perc_mns<-ggplot(data = bf_comb, aes(x=GW, y=PVpercsum_mn, group=Fishery, color = Fishery, fill = Fishery)) + # Nominal expected impact across targets
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = PVlperc, ymax = PVuperc),linetype=2, alpha=0.3) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 11, y = min(PVpercsum_mn) - diff(range(PVpercsum_mn)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_label(aes(x = 55, y = min(PVpercsum_mn) - diff(range(PVpercsum_mn)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean Percent", colour = "Fishery") +
  xlim(c(0,60))

# Relabeling plots
bf_nom_md<-bf_nom_md + labs(y = "Median present value exposed ($Mil)")
bf_nom_mdCA<-bf_nom_mdCA + labs(y = "Median present value exposed ($Mil)")
bf_nom_mdOR<-bf_nom_mdOR + labs(y = "Median present value exposed ($Mil)")
bf_nom_mds<-bf_nom_mds + labs(y = "Median present value exposed ($Mil)")
bf_nom_mdWA<-bf_nom_mdWA + labs(y = "Median present value exposed ($Mil)")
bf_nom_mn<-bf_nom_mn + labs(y = "Mean present value exposed ($Mil)") + xlim(c(0,65))
bf_nom_mnCA<-bf_nom_mnCA + labs(y = "Mean present value exposed ($Mil)")
bf_nom_mnOR<-bf_nom_mnOR + labs(y = "Mean present value exposed ($Mil)")
bf_nom_mns<-bf_nom_mns + labs(y = "Mean present value exposed ($Mil)") + xlim(c(0,65))
bf_nom_mnWA<-bf_nom_mnWA + labs(y = "Mean present value exposed ($Mil)")
bf_perc_md<-bf_perc_md + labs(y = "Median present value exposed (% of fishery value)")
bf_perc_mdCA<-bf_perc_mdCA + labs(y = "Median present value exposed (% of fishery value)")
bf_perc_mdOR<-bf_perc_mdOR + labs(y = "Median present value exposed (% of fishery value)")
bf_perc_mds<-bf_perc_mds + labs(y = "Median present value exposed (% of fishery value)")
bf_perc_mdWA<-bf_perc_mdWA + labs(y = "Median present value exposed (% of fishery value)")
bf_perc_mn<-bf_perc_mn + labs(y = "Mean present value exposed (% of fishery value)") + xlim(c(0,65))
bf_perc_mnCA<-bf_perc_mnCA + labs(y = "Mean present value exposed (% of fishery value)")
bf_perc_mnOR<-bf_perc_mnOR + labs(y = "Mean present value exposed (% of fishery value)")
bf_perc_mns<-bf_perc_mns + labs(y = "Mean present value exposed (% of fishery value)") + xlim(c(0,65))
bf_perc_mnWA<-bf_perc_mnWA + labs(y = "Mean present value exposed (% of fishery value)")

ci_bf_nom_mn<-ci_bf_nom_mn + labs(y = "Mean present value exposed ($Mil)") + xlim(c(0,65))
ci_bf_nom_mnCA<-ci_bf_nom_mnCA + labs(y = "Mean present value exposed ($Mil)")
ci_bf_nom_mnOR<-ci_bf_nom_mnOR + labs(y = "Mean present value exposed ($Mil)")
ci_bf_nom_mns<-ci_bf_nom_mns + labs(y = "Mean present value exposed ($Mil)") + xlim(c(0,65))
ci_bf_nom_mnWA<-ci_bf_nom_mnWA + labs(y = "Mean present value exposed ($Mil)")
ci_bf_perc_mn<-ci_bf_perc_mn + labs(y = "Mean present value exposed (% of fishery value)") + xlim(c(0,65))
ci_bf_perc_mnCA<-ci_bf_perc_mnCA + labs(y = "Mean present value exposed (% of fishery value)")
ci_bf_perc_mnOR<-ci_bf_perc_mnOR + labs(y = "Mean present value exposed (% of fishery value)")
ci_bf_perc_mns<-ci_bf_perc_mns + labs(y = "Mean present value exposed (% of fishery value)") + xlim(c(0,65))
ci_bf_perc_mnWA<-ci_bf_perc_mnWA + labs(y = "Mean present value exposed (% of fishery value)")

# Grouped plots
t_nom_r_mn<-bf_nom_mn + bf_nom_mns + plot_layout(axis_titles = "collect") 
t_nom_s_mn<-bf_nom_mnCA + bf_nom_mnOR + bf_nom_mnWA + plot_layout(axis_titles = "collect")
t_nom_r_mn / t_nom_s_mn + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")

t_nom_r_md<-bf_nom_md + bf_nom_mds + plot_layout(axis_titles = "collect") 
t_nom_s_md<-bf_nom_mdCA + bf_nom_mdOR + bf_nom_mdWA + plot_layout(axis_titles = "collect")
t_nom_r_md / t_nom_s_md + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")

t_perc_r_md<-bf_perc_md + bf_perc_mds + plot_layout(axis_titles = "collect") 
t_perc_s_md<-bf_perc_mdCA + bf_perc_mdOR + bf_perc_mdWA + plot_layout(axis_titles = "collect")
t_perc_r_md / t_perc_s_md + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")

t_perc_r_mn<-bf_perc_mn + bf_perc_mns + plot_layout(axis_titles = "collect") 
t_perc_s_mn<-bf_perc_mnCA + bf_perc_mnOR + bf_perc_mnWA + plot_layout(axis_titles = "collect")
t_perc_r_mn / t_perc_s_mn + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")

ci_t_nom_r_mn<-ci_bf_nom_mn + ci_bf_nom_mns + plot_layout(axis_titles = "collect") 
ci_t_nom_s_mn<-ci_bf_nom_mnCA + ci_bf_nom_mnOR + ci_bf_nom_mnWA + plot_layout(axis_titles = "collect")
ci_t_nom_r_mn / ci_t_nom_s_mn + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")

ci_t_perc_r_mn<-ci_bf_perc_mn + ci_bf_perc_mns + plot_layout(axis_titles = "collect") 
ci_t_perc_s_mn<-ci_bf_perc_mnCA + ci_bf_perc_mnOR + ci_bf_perc_mnWA + plot_layout(axis_titles = "collect")
ci_t_perc_r_mn / ci_t_perc_s_mn + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")

# Combined groups plot

# layout<-"AAABBB
#          CCDDEE
#          FFFGGG
#          HHIIJJ"
# 
# layout <- c(
#   area(t = 1, l = 1, b = 1, r = 3),
#   area(t = 1, l = 4, b = 1, r = 6),
#   area(t = 2, l = 1, b = 2, r = 2),
#   area(t = 2, l = 3, b = 2, r = 4),
#   area(t = 2, l = 5, b = 2, r = 6),
#   area(t = 3, l = 1, b = 3, r = 3),
#   area(t = 3, l = 4, b = 3, r = 6),
#   area(t = 4, l = 1, b = 4, r = 2),
#   area(t = 4, l = 3, b = 4, r = 4),
#   area(t = 4, l = 5, b = 4, r = 6)
# )

layout<-"ABCDE
         "
big_nom_mn<-bf_nom_mn + bf_nom_mns + bf_nom_mnCA + bf_nom_mnOR + bf_nom_mnWA + plot_layout(axis_titles = "collect", design = layout)
big_perc_mn<-bf_perc_mn + bf_perc_mns + bf_perc_mnCA + bf_perc_mnOR + bf_perc_mnWA + plot_layout(axis_titles = "collect", design = layout)
(big_nom_mn / big_perc_mn) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')

p4nom<-ggplot(data.frame(l = "Mean present value exposed ($Mil)", x = 1, y = 1)) + # Collecting axes titles across rows https://stackoverflow.com/questions/65291723/merging-two-y-axes-titles-in-patchwork
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

p4perc<-ggplot(data.frame(l = "Mean present value exposed (% of fishery value)", x = 1, y = 1)) + # Collecting axes titles across rows https://stackoverflow.com/questions/65291723/merging-two-y-axes-titles-in-patchwork
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

p5<-ggplot(data.frame(l = bf_nom_mn$labels$x, x = 1, y = 1)) + # Collecting axes titles across rows https://stackoverflow.com/questions/65291723/merging-two-y-axes-titles-in-patchwork
  geom_text(aes(x, y, label = l), angle = 0) + 
  theme_void() +
  coord_cartesian(clip = "off")

bf_nom_mn$labels$y<-bf_nom_mns$labels$y<-bf_nom_mnCA$labels$y<-bf_nom_mnOR$labels$y<-bf_nom_mnWA$labels$y<-""
bf_nom_mn$labels$x<-bf_nom_mns$labels$x<-bf_nom_mnCA$labels$x<-bf_nom_mnOR$labels$x<-bf_nom_mnWA$labels$x<-""
bf_perc_mn$labels$y<-bf_perc_mns$labels$y<-bf_perc_mnCA$labels$y<-bf_perc_mnOR$labels$y<-bf_perc_mnWA$labels$y<-""
bf_perc_mn$labels$x<-bf_perc_mns$labels$x<-bf_perc_mnCA$labels$x<-bf_perc_mnOR$labels$x<-bf_perc_mnWA$labels$x<-""

lay<-"ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ABBBBBBBBBCCCCCCCCC
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      ADDDDDDEEEEEEFFFFFF
      GGGGGGGGGGGGGGGGGGG
      "
big_nom_mn<-(p4nom + bf_nom_mn + bf_nom_mns + bf_nom_mnCA + bf_nom_mnOR + bf_nom_mnWA + p5) + plot_layout(design = lay)
big_perc_mn<-(p4perc + bf_perc_mn + bf_perc_mns + bf_perc_mnCA + bf_perc_mnOR + bf_perc_mnWA + p5) + plot_layout(design = lay)
(big_nom_mn / big_perc_mn) + plot_layout(guides = "collect") + plot_annotation(tag_levels =  list(c("","A","B","C","D","E","","","F","G","H","I","J",""))) # 1200*1400 export

big_nom_mn<-((bf_nom_mn + bf_nom_mns) / (bf_nom_mnCA + bf_nom_mnOR + bf_nom_mnWA)) + plot_layout(axis_titles = "collect", guides = "collect")
big_perc_mn<-(bf_perc_mn + bf_perc_mns + bf_perc_mnCA + bf_perc_mnOR + bf_perc_mnWA) + plot_layout(axis_titles = "collect", design = lay)

## Total regional exposure in 2045 when optimizing regionally or by state 
bf_comb %>% filter(GW == 59.4) %>% summarise(sumPV = sum(PVnomsum_mn)/1000000)
biggerfishR %>% filter(n==62) %>% summarise(sumPV = sum(`Mean PV`)/1000000)


# Unused code testing alternate optimization approaches -------------------
# ## Non-domiNULL# ## Non-dominated Sorting Genetic Algorithm II using a repair function and index based drawing
# 
# #df<-df[sample(nrow(df), replace = FALSE, size = nrow(df)/((300/8.185))), ] # Testing with smaller dataframes D=9 for region (technically 8.8888; 8000/900) when cells are 300km2
# 
# # ggplot() + # Trying to remove clearly dominated locations
# #   geom_point(data = df2, aes(x = Total_fish_USD, y = LCOE_MWh, color = "Red"), size = 2) +
# #   theme_minimal()
# # 
# # df2<-df[df$LCOE_MWh<200&df$Total_fish_USD<10,]
# 
# optimized_repair_solution <- function(solution, n_select, n_items) {
#   # Ensure the solution does not exceed the desired length, and remove duplicates
#   unique_solution <- unique(round(solution),0)
#   
#   # If the unique solution has more items than needed, randomly select n_select items
#   if(length(unique_solution) > n_select) {
#     return(sample(unique_solution, n_select))
#   }
#   
#   # If the unique solution has fewer items than needed, fill in with missing items
#   if(length(unique_solution) < n_select) {
#     # Pre-calculate all missing items once
#     all_items <- 1:n_items
#     missing_items <- all_items[!all_items %in% unique_solution]
#     # Randomly select enough missing items to reach n_select
#     to_add <- sample(missing_items, n_select - length(unique_solution))
#     unique_solution <- c(unique_solution, to_add)
#   }
#   
#   # Return the repaired solution, ensuring it's sorted for consistency
#   return(sort(unique_solution))
# }
# 
# main.goal1<-function(x){    # x - a vector of row indices
#   x_repaired <- optimized_repair_solution(x, D, nrow(df))
#   return(sum(df[x_repaired, 4]))}
# 
# main.goal2<-function(x){    # x - a vector of row indices
#   x_repaired <- optimized_repair_solution(x, D, nrow(df))
#   return(mean(df[x_repaired, 6]))}
# 
# eval<-function(x){
#   return(c(main.goal1(x),main.goal2(x)))} # objective function
# 
# D<-9 # Number of values in the function must equal D
# 
# system.time(G<-nsga2(fn=eval,
#          idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn 
#          odim=2, # Output dimensions
#          lower.bounds=rep(0,nrow(df)),
#          upper.bounds=rep(as.numeric(nrow(df)),D),
#          popsize=200,generations=1000, cprob = 0.7, cdist = 5,
#          mprob = 0.2, mdist = 10))
# plot(G)
# plot(G$par)
# 
# # Testing if optimal at region scale
# llcoe<-lowest_n_df <- df[order(df$LCOE_MWh)[1:D], ]
# sum(llcoe$Total_fish_USD)
# mean(llcoe$LCOE_MWh)
# 
# llcoe_pt<-as.data.frame(cbind(sum(llcoe$Total_fish_USD),mean(llcoe$LCOE_MWh)))
# names(llcoe_pt) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)") # Rename columns
# 
# lfpv<-lowest_n_df <- df[order(df$Total_fish_USD)[1:D], ]
# sum(lfpv$Total_fish_USD)
# mean(lfpv$LCOE_MWh)
# 
# lfpv_pt<-as.data.frame(cbind(sum(lfpv$Total_fish_USD),mean(lfpv$LCOE_MWh)))
# names(lfpv_pt) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)") # Rename columns
# 
# # Plotting
# pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
# names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
# pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
# pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]
# 
#   
# a<-ggplot() +
#   geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
#   geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
#   #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
#   #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
#   #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
#   theme_minimal() +
#   #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
#   theme(legend.title = element_blank(),legend.position="none") 

## Below creates a matrix of unique combinations of values of length m. Works for very small vectors, quickly becomes too big to store in R with our vector length. Going to need to optimize or sample.
#dft<-RcppAlgos::comboGeneral(as.vector(df$`Wind farm grid ID`),m = D, repetition = FALSE) # https://stackoverflow.com/questions/65007670/create-unique-possible-combinations-from-the-elements-in-vector-in-r

# x<-df
# it<-1000
# n<-326

# ## Sampling based approach
# pareto<-function(x,n,it){
#   # x - Dataframe
#   # n - Subset size (# of cells needed to meet wind energy goal)
#   # it - Number of sampling iterations
#   
#   # Storage for all evaluated solutions
#   all_solutions <- data.frame(sum_v1 = numeric(), avg_v2 = numeric(), CA = numeric(), WA = numeric(), OR = numeric(), sum_MWhr_yr = numeric())
#   
#   for (i in 1:it) {
#     # Randomly select n items
#     selected_indices <- sample(nrow(x), n)
#     
#     # Calculate objectives
#     sum_v1 <- sum(x$Total_fish_USD[selected_indices])  # Objective 1: Minimize sum of v1
#     avg_v2 <- mean(x$LCOE_MWh[selected_indices])  # Objective 2: Minimize average of v2
#     
#     # N in each state and aggregate energy
#     a<-x %>% slice(selected_indices)
#     e<-sum(a$MWhyraw) # Aggregate energy
#     a<-a %>% 
#       group_by(State) %>% 
#       summarise(count = n(),.groups = "drop") %>% 
#       pivot_wider(names_from = State,values_from = count)
# 
#     # Store each solution's objectives
#     all_solutions <- rbind(all_solutions, list(sum_v1 = sum_v1, avg_v2 = avg_v2, CA = a$CA, WA = a$WA, OR = a$OR, sum_MWhr_yr = e))
#   }
#   
#   # Determine Pareto-optimality for each solution
#   is_pareto_optimal <- function(index, solutions) {
#     for (j in 1:nrow(solutions)) {
#       if (j != index && solutions[j, "sum_v1"] <= solutions[index, "sum_v1"] && solutions[j, "avg_v2"] <= solutions[index, "avg_v2"]) {
#         return(FALSE)  # Solution is dominated if another solution is better (lower) in both objectives
#       }
#     }
#     return(TRUE)  # Solution is Pareto-optimal if no other solution is better (lower) in both objectives
#   }
#   
#   pareto_optimal_indices <- sapply(1:nrow(all_solutions), function(i) is_pareto_optimal(i, all_solutions))
#   
#   # Extract Pareto-optimal solutions
#   pareto_solutions <- all_solutions[pareto_optimal_indices, ]
#   pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$sum_v1), ]
#   
#   # Adding an identifier for plotting
#   all_solutions$Type <- 'Dominated'
#   pareto_solutions_ordered$Type <- 'Pareto-optimal'
#   
#   # Combine all solutions for plotting
#   plot_data <- rbind(
#     all_solutions,
#     pareto_solutions_ordered)
#   
#   a<-ggplot() +
#     geom_point(data = plot_data %>% filter(Type=="Pareto-optimal"), aes(x = sum_v1, y = avg_v2, color = "Red"), size = 2) +
#     geom_line(data = pareto_solutions_ordered, aes(x = sum_v1, y = avg_v2), color = "blue", linewidth = 1) +
#     #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
#     theme_minimal() +
#     labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
#     theme(legend.title = element_blank(),legend.position="none") 
#   
#   return(list(a,plot_data))
# }
# 
# paretos<-function(x,n,it){ # Same as function above, but doesn't track outcomes by state
#   # x - Dataframe
#   # n - Subset size (# of cells needed to meet wind energy goal)
#   # it - Number of sampling iterations
#   
#   # Storage for all evaluated solutions
#   all_solutions <- data.frame(sum_v1 = numeric(), avg_v2 = numeric(), sum_MWhr_yr = numeric())
#   
#   for (i in 1:it) {
#     # Randomly select n items
#     selected_indices <- sample(nrow(x), n)
#     
#     # Calculate objectives
#     sum_v1 <- sum(x$Total_fish_USD[selected_indices])  # Objective 1: Minimize sum of v1
#     avg_v2 <- mean(x$LCOE_MWh[selected_indices])  # Objective 2: Minimize average of v2
#     
#     # Aggregate energy 
#     a<-x %>% slice(selected_indices)
#     e<-sum(a$MWhyraw) # Aggregate energy
#     
#     # Store each solution's objectives
#     all_solutions <- rbind(all_solutions, list(sum_v1 = sum_v1, avg_v2 = avg_v2, sum_MWhr_yr = e))
#   }
#   
#   # Determine Pareto-optimality for each solution
#   is_pareto_optimal <- function(index, solutions) {
#     for (j in 1:nrow(solutions)) {
#       if (j != index && solutions[j, "sum_v1"] <= solutions[index, "sum_v1"] && solutions[j, "avg_v2"] <= solutions[index, "avg_v2"]) {
#         return(FALSE)  # Solution is dominated if another solution is better (lower) in both objectives
#       }
#     }
#     return(TRUE)  # Solution is Pareto-optimal if no other solution is better (lower) in both objectives
#   }
#   
#   pareto_optimal_indices <- sapply(1:nrow(all_solutions), function(i) is_pareto_optimal(i, all_solutions))
#   
#   # Extract Pareto-optimal solutions
#   pareto_solutions <- all_solutions[pareto_optimal_indices, ]
#   pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$sum_v1), ]
#   
#   # Adding an identifier for plotting
#   all_solutions$Type <- 'Dominated'
#   pareto_solutions_ordered$Type <- 'Pareto-optimal'
#   
#   # Combine all solutions for plotting
#   plot_data <- rbind(
#     all_solutions,
#     pareto_solutions_ordered)
#   
#   a<-ggplot() +
#     geom_point(data = plot_data %>% filter(Type=="Pareto-optimal"), aes(x = sum_v1, y = avg_v2, color = "Red"), size = 2) +
#     geom_line(data = pareto_solutions_ordered, aes(x = sum_v1, y = avg_v2), color = "blue", linewidth = 1) +
#     #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
#     theme_minimal() +
#     labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
#     theme(legend.title = element_blank(),legend.position="none") 
#   
#   return(list(a,plot_data))
# }
# 
# # Sampling
# iterations<-1000
# system.time(region<-pareto(x = df,n = 326, it = iterations)) # 9 hrs @ 1m iterations
# system.time(ca<-paretos(x = df %>% filter(State == "CA"),n = 204, it = iterations)) # 3.5 hrs @ 1m
# system.time(or<-paretos(x = df %>% filter(State == "OR"),n = 122, it = iterations)) # 3.5 hrs @ 1m
# system.time(wa<-paretos(x = df %>% filter(State == "WA"),n = 122, it = iterations)) # 3.5 hrs @ 1m
# 
# (region[[1]] + ca[[1]]) / (or[[1]] + wa[[1]]) + # Unified plot
#   plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 12))
# 
# region[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal") # Number of pareto optimal sites in each state across scenarios
# 
# region[[2]]$lcoe_x_MWh<-region[[2]]$avg_v2*region[[2]]$sum_MWhr_yr # Cost of production 
# ca[[2]]$lcoe_x_MWh<-ca[[2]]$avg_v2*ca[[2]]$sum_MWhr_yr
# wa[[2]]$lcoe_x_MWh<-wa[[2]]$avg_v2*wa[[2]]$sum_MWhr_yr
# or[[2]]$lcoe_x_MWh<-or[[2]]$avg_v2*or[[2]]$sum_MWhr_yr
# 
# region[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal") # Number of pareto optimal sites in each state across scenarios
# 
# ca[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal") # Table of pareto optimal outcomes
# or[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal")
# 
# df %>% group_by(State) %>% # Mean LCOE by state
#   summarise(meanLCOE = mean(LCOE_MWh))
# 
# ggplot(df, aes(LCOE_MWh, colour = State)) + # CDF of LCOE by state
#   stat_ecdf(linewidth = 1) +
#   labs(y = "", x = "Mean LCOE ($/MWh)") +
#   theme_minimal()
# 
# # ggplot(df, aes(LCOE_MWh, weight = weight, color = State)) + # Histogram by state to account for differing sample sizes
# #   #geom_density() +
# #   stat_bin(geom = "line", bins = 500, position = position_identity(),linewidth = 1) +
# #   scale_fill_brewer(palette = "Pastel1") +
# #   theme_minimal()
# 
# 
# ## GPareto - probably need to adapt the function more for our context
# 
# f<-function(x){
#   print(nrow(x)) # Debug line to print the number of rows in x
#   if (!is.data.frame(x)) {
#     stop("The input x is not a dataframe.")
#   }
#   if (nrow(x) < 3) {
#     stop("The dataframe x must have at least 3 rows.")
#   }
#   s<-x[sample(nrow(x), replace = FALSE, size = 3), ]
#   sx<-sum(s$Total_fish_USD)
#   sy<-mean(s$LCOE_kWh)
#   
#   print(s)
#   print(sx)
#   print(sy)
#   
#   return(cbind(sx,sy))
# }
# 
# f2<-function(x){
#   s<-as.matrix(x)
#   s<-s[sample(nrow(s), replace = FALSE, size = 3), , drop = FALSE] ## See this https://stackoverflow.com/questions/50843515/sampling-columns-in-a-matrix-within-r about dropping
#   sx<-sum(s[,1])
#   sy<-mean(s[,2], na.rm = TRUE)
#   return(cbind(sx,sy))
# }
# 
# f3<-function(x){
#   s<-df3[sample(nrow(df3), replace = FALSE, size = 3),]
#   sx<-sum(s[,1])
#   sy<-mean(s[,2], na.rm = TRUE)
#   return(cbind(sx,sy))
# }
# 
# f2(df3)
# 
# df3<-as.matrix(df3)
# 
# res<-easyGParetoptim(fn = f2, budget = 5, lower = c(0,0), upper = c(1,1))
# 
# 
# ####
# P1 <-
#   function(x){
#     if(is.null(dim(x))){
#       x <- matrix(x, nrow = 1) 
#     }
#     b1<-15*x[,1]-5
#     b2<-15*x[,2]
#     return(cbind((b2-5.1*(b1/(2*pi))^2+5/pi*b1-6)^2 +10*((1-1/(8*pi))*cos(b1)+1),
#                  -sqrt((10.5-b1)*(b1+5.5)*(b2+0.5)) - 1/30*(b2 -5.1*(b1/(2*pi))^2-6)^2 - 1/3*((1-1/(8*pi))*cos(b1)+1)
#     ) 
#     )
#   }
# 
# d <- 2; ninit <- 10; fun <- P1
# design <- DiceDesign::lhsDesign(ninit, d, seed = 42)$design
# class(design)
# 
# res<-easyGParetoptim(fn = P1, budget = 50, lower = c(0,0), upper = c(1,1))
# plotGPareto(res)
# 
# #####
# MOP2 <- function(x)
# {
#   xmod <- x*4 - 2
#   if (is.null(nrow(x)))
#   { 
#     n <- length(xmod)
#     y1 <- 1 - exp(-sum((xmod - 1/sqrt(n))^2) )
#     y2 <- 1 - exp(-sum((xmod + 1/sqrt(n))^2) )
#     Y <- matrix(c(y1,y2),1,2)
#   } else
#   {
#     n <- ncol(xmod)
#     y1 <- 1 - exp(-rowSums((xmod - 1/sqrt(n))^2) )
#     y2 <- 1 - exp(-rowSums((xmod + 1/sqrt(n))^2) )
#     Y <- cbind(y1,y2)
#   }
#   
#   return(Y)
# }
# 
# xmod <- a*4 - 2
# 
# design.init <- matrix(seq(0, 1, length.out = 6), ncol = 1)
# response.init <- MOP2(design.init)
# MOP2(design.init)
# 
# mf1 <- km(~1, design = design.init, response = response.init[, 1])
# mf2 <- km(~1, design = design.init, response = response.init[, 2])
# model <- list(mf1, mf2)
# 
# res <- GParetoptim(model = model, fn = MOP2, crit = "EHI", nsteps = 25,lower = 0, upper = 1, critcontrol = list(refPoint = c(2, 2)))
# plotGPareto(res)
# res<-easyGParetoptim(fn = MOP2, budget = 100, lower = c(0,0), upper = c(1,1))
# plotGPareto(res)
# 
# 
# ###############
# 
# DTLZ2 <- function(x, nobj = 3){
#   if(is.null(dim(x))){
#     x <- matrix(x, 1) 
#   }
#   n <- ncol(x)
#   
#   y <- matrix(x[,1:(nobj-1)], nrow(x))
#   z <- matrix(x[,nobj:n], nrow(x))
#   
#   g <- rowSums((z-0.5)^2)
#   
#   #   tmp <- c(rev(cumprod(cos(y * pi/2))), 1)
#   #   tmp2 <- c(1, rev(sin(y * pi/2)))
#   tmp <- t(apply(cos(y * pi/2), 1, cumprod))
#   tmp <- cbind(t(apply(tmp, 1, rev)), 1)
#   
#   tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
#   
#   f <- tmp * tmp2 * (1 + g)
#   
# }
# 
# y <- matrix(data = x[,1:(3-1)], nrow(x))
# z <- matrix(x[,3:n], nrow(x))
# 
# a<-rbind(rep(1,5),seq(1,5,1))
# 
# 
# 
# DTLZ2
# 
# ZDT2(a)
# 
# ?ZDT3()
# 
#     
