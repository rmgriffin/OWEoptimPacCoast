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

## Loading data (all locations at 8.185 km2 per cell)
# df<-read_csv("OWEP LCOE NPV & fishing PV data V2 DO NOT DISTRIBUTE.csv")
# df<-df %>%
#   dplyr::select(OWEP_ID,Total_fish_USD,LCOE_kWh,State,MWhyraw,Gridcell_area_m2)
df<-read_csv("OWEP output & fishing PV data V4 DO NOT DISTRIBUTE.csv")
df<-df %>%
  dplyr::select(`Wind farm grid ID`,`Wind farm state`,`OWEP grid cells (n)`,Total_fish_USD,`Weighted mean LCOE`,`Total MWhyraw`,`Wind farm area (m2)`)

#sapply(df, function(y) sum(length(which(is.na(y))))) # Lots of locations with fishing not suitable for wind presumably, lots of NAs for LCOE
df<-df %>% drop_na(`Weighted mean LCOE`)
#sum(df$Gridcell_area_m2)/1000000 # Total area over which to optimize version 2

df<-df[df$`OWEP grid cells (n)`>34,] # Choosing 34 as a cutoff for what could be considered full cells, could amend
#sum(df$`Wind farm area (m2)`)/1000000 # Total area over which to optimize version 4

df$Total_fish_USD<-df$Total_fish_USD/1000000
df$LCOE_MWh<-df$`Weighted mean LCOE`*1000
df$`Weighted mean LCOE`<-NULL
df$`Wind farm area (m2)`<-NULL # This isn't a correct calculation, this is actually cell size
df<-as.data.frame(df)
df2<-df

## State targets. CA is 5GW by 2030 and 25GW by 2045, OR is 3GW by 2030. Could simulate WA at 3GW and 2045 by multiplying all by 5 (Region 63, CA 28, OR 17, WA 17)
# Assuming 3MW per km2 and each cell is 8.185km2, CA is 1667km2 (204 cells) by 2030, OR is 1000km2 (122 cells) by 2030, and together is 2667km2 (326 cells). 
# Assuming 3MW per km2 and each cell is 295km2 (36*8.185km2 cells), CA is 1667km2 (6 cells) by 2030, OR is 1000km2 (4 cells) by 2030, assuming WA is 1000km2 (4 cells) by 2030, and together is 3667km2 (13 cells) 

## Non-dominated Sorting Genetic Algorithm II using binary based drawing and a repair function to set constraint (penalties and the constraint function don't seem to work)

transform_vector_to_exact_ones <- function(v, num_ones) { # Sets constraint by forcing the NSGA2 draws to generate exactly n ones and (D - n) zeros
  # Length of the input vector
  n <- length(v)
  
  # Create a vector of zeros
  transformed_v <- rep(0, n)
  
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

# Region ---

D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-13 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=500,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

#plot(G)

# Identifying pareto sets after optimization
pareto_indices <- which(G$pareto.optimal == TRUE)
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions$Pareto<-NULL
pareto_solutions<-unique(pareto_solutions)

transform_row_to_exact_ones <- function(row, num_ones) { # Same as vector version above, but for a matrix
  n <- length(row)
  transformed_row <- rep(0, n)
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
selected_wfgridID<-as.data.frame(selected_rows) %>%
  rowwise() %>%
  mutate(across(everything(), ~df$`Wind farm grid ID`[.x])) %>%
  ungroup()
agree<-as.data.frame(table(as.matrix(selected_wfgridID)))
agree$Var1<-as.numeric(levels(agree$Var1))[agree$Var1] # https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
agree$prob<-agree$Freq/sum(agree$Freq)

# Fisheries exposure at pareto selected sites (sampling approach, sample n from agreement matrix weighted by how often they are selected and calculate an expected value)
dft<-read_csv("OWEP output & fishing PV data V4 DO NOT DISTRIBUTE.csv")
dft<-dft %>% dplyr::select(`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD)
dft<-as.data.frame(dft)

m<-1000
bigfish<-data.frame("Dungeness_USD"=rep(NA,m),"At-sea_hake_USD"=rep(NA,m),"Shore_hake_USD"=rep(NA,m),"Market_squid_USD"=rep(NA,m),"Pink_shrimp_USD"=rep(NA,m),"Albacore_USD"=rep(NA,m),"Chinook_USD"=rep(NA,m),"Sablefish_USD"=rep(NA,m),"Spiny_lobster_USD"=rep(NA,m),"Sites"=rep(NA,m))
for(j in 1:m){ # Sums species PV for a drawn sample of size n and iterates across samples while tracking draws in "Sites" column
  s_agree<-sample(agree$Var1,n,replace = FALSE,prob = agree$prob)
  fishroll<-data.frame("Dungeness_USD"=rep(NA,n),"At-sea_hake_USD"=rep(NA,n),"Shore_hake_USD"=rep(NA,n),"Market_squid_USD"=rep(NA,n),"Pink_shrimp_USD"=rep(NA,n),"Albacore_USD"=rep(NA,n),"Chinook_USD"=rep(NA,n),"Sablefish_USD"=rep(NA,n),"Spiny_lobster_USD"=rep(NA,n))
  for(i in 1:n){ # Identifies species PV from a given weighted sample of size n drawn from agreement vector
    fishroll[i,]<-dft[dft$`Wind farm grid ID`==s_agree[i],2:10]
  }
  fishroll<-as.data.frame(t(colSums(fishroll))) #paste(s_agree, collapse = " ")
  fishroll<-cbind(fishroll,paste(s_agree, collapse = " "))
  names(fishroll)[names(fishroll) == 'paste(s_agree, collapse = " ")']<-"Sites"
  bigfish[j,]<-fishroll
}

bigfish<-bigfish %>% pivot_longer(cols = !Sites,names_to = "Fishery",values_to = "PV")
bigfish$Fishery<-gsub("_USD$", "",bigfish$Fishery)
bigfish$Fishery<-gsub("_", " ",bigfish$Fishery)
bigfish$Fishery<-gsub("\\.", "-",bigfish$Fishery)

sorted_fish<-bigfish %>% distinct(Fishery) %>% 
  arrange(desc(`Fishery`)) %>% unlist()

bigfish$Fishery<-factor(bigfish$Fishery, levels = sorted_fish)

bigfish30<-ggplot() + 
  geom_boxplot(data = bigfish, aes(y = Fishery, x = PV/1000000),outlier.shape = NA) +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)", color = "Fishery") +
  labs(y="",x = "Fishing PV ($Mil)") +
  coord_cartesian(xlim = c(0,7)) +
  theme_minimal()

fishsum<-dft %>% pivot_longer(cols = !`Wind farm grid ID`,names_to = "Fishery",values_to = "PV") %>% 
  group_by(Fishery) %>% 
  summarise(fishsum = sum(PV))

fishsum$Fishery<-gsub("_USD$", "",fishsum$Fishery)
fishsum$Fishery<-gsub("_", " ",fishsum$Fishery)
fishsum$Fishery<-gsub("\\.", "-",fishsum$Fishery)

bigfish<-merge(bigfish,fishsum,by="Fishery")
bigfish$PVperc<-bigfish$PV/bigfish$fishsum

bigfish30perc<-ggplot() + 
  geom_boxplot(data = bigfish, aes(y = Fishery, x = PVperc*100),outlier.shape = NA) +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)", color = "Fishery") +
  labs(y="",x = "Pecentage of total Fishing PV (%)") +
  coord_cartesian(xlim = c(0,1)) +
  theme_minimal()

rm(bigfish,dft,agree,fishroll,idim_par,selected_wfgridID,selected_rows,i,j,m,pareto_indices,s_agree,sorted_fish,fishsum)

# Testing if optimal at region scale
llcoe<-df[order(df$LCOE_MWh)[1:n], ]
sum(llcoe$Total_fish_USD)
mean(llcoe$LCOE_MWh)

llcoe_pt<-as.data.frame(cbind(sum(llcoe$Total_fish_USD),mean(llcoe$LCOE_MWh)))
names(llcoe_pt) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)") # Rename columns

lfpv<-df[order(df$Total_fish_USD)[1:n], ]
sum(lfpv$Total_fish_USD)
mean(lfpv$LCOE_MWh)

lfpv_pt<-as.data.frame(cbind(sum(lfpv$Total_fish_USD),mean(lfpv$LCOE_MWh)))
names(lfpv_pt) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)") # Rename columns

# Plotting 
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

region<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none") +
  coord_cartesian(xlim = c(0, 80))

# CA ---
df<-df2[df2$`Wind farm state`=="CA",]
D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-6 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=200,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

# Plotting
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

CA<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none") +
  coord_cartesian(xlim = c(0, 80))

# OR ---
df<-df2[df2$`Wind farm state`=="OR",]
D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-4 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=200,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

# Plotting
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

OR<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none")

# WA ---
df<-df2[df2$`Wind farm state`=="WA",]
D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-4 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=200,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

# Plotting
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

WA<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none")

# Region 2045 ---
df<-df2
D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-63 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=200,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

# Identifying pareto sets after optimization
pareto_indices <- which(G$pareto.optimal == TRUE)

transform_row_to_exact_ones <- function(row, num_ones) { # Same as vector version above, but for a matrix
  n <- length(row)
  transformed_row <- rep(0, n)
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
selected_rows<-unique(selected_rows) # NSGA2 can produce duplicate pareto sets
selected_wfgridID<-as.data.frame(selected_rows) %>%
  rowwise() %>%
  mutate(across(everything(), ~df$`Wind farm grid ID`[.x])) %>%
  ungroup()
agree<-as.data.frame(table(as.matrix(selected_wfgridID)))
agree$Var1<-as.numeric(levels(agree$Var1))[agree$Var1] # https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
agree$prob<-agree$Freq/sum(agree$Freq)

# Fisheries exposure at pareto selected sites (sampling approach, sample n from agreement matrix weighted by how often they are selected and calculate an expected value)
dft<-read_csv("OWEP output & fishing PV data V4 DO NOT DISTRIBUTE.csv")
dft<-dft %>% dplyr::select(`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD)
dft<-as.data.frame(dft)

m<-10000 # Higher than for 2030 because there are more potential sites in 2045 (119 sites are observed across pareto sets)
bigfish<-data.frame("Dungeness_USD"=rep(NA,m),"At-sea_hake_USD"=rep(NA,m),"Shore_hake_USD"=rep(NA,m),"Market_squid_USD"=rep(NA,m),"Pink_shrimp_USD"=rep(NA,m),"Albacore_USD"=rep(NA,m),"Chinook_USD"=rep(NA,m),"Sablefish_USD"=rep(NA,m),"Spiny_lobster_USD"=rep(NA,m),"Sites"=rep(NA,m))
for(j in 1:m){ # Sums species PV for a drawn sample of size n and iterates across samples while tracking draws in "Sites" column
  s_agree<-sample(agree$Var1,n,replace = FALSE,prob = agree$prob)
  fishroll<-data.frame("Dungeness_USD"=rep(NA,n),"At-sea_hake_USD"=rep(NA,n),"Shore_hake_USD"=rep(NA,n),"Market_squid_USD"=rep(NA,n),"Pink_shrimp_USD"=rep(NA,n),"Albacore_USD"=rep(NA,n),"Chinook_USD"=rep(NA,n),"Sablefish_USD"=rep(NA,n),"Spiny_lobster_USD"=rep(NA,n))
  for(i in 1:n){ # Identifies species PV from a given weighted sample of size n drawn from agreement vector
    fishroll[i,]<-dft[dft$`Wind farm grid ID`==s_agree[i],2:10]
  }
  fishroll<-as.data.frame(t(colSums(fishroll))) #paste(s_agree, collapse = " ")
  fishroll<-cbind(fishroll,paste(s_agree, collapse = " "))
  names(fishroll)[names(fishroll) == 'paste(s_agree, collapse = " ")']<-"Sites"
  bigfish[j,]<-fishroll
}

bigfish<-bigfish %>% pivot_longer(cols = !Sites,names_to = "Fishery",values_to = "PV")
bigfish$Fishery<-gsub("_USD$", "",bigfish$Fishery)
bigfish$Fishery<-gsub("_", " ",bigfish$Fishery)
bigfish$Fishery<-gsub("\\.", "-",bigfish$Fishery)

sorted_fish<-bigfish %>% distinct(Fishery) %>% 
  arrange(desc(`Fishery`)) %>% unlist()

bigfish$Fishery<-factor(bigfish$Fishery, levels = sorted_fish)

bigfish45<-ggplot() + 
  geom_boxplot(data = bigfish, aes(y = Fishery, x = PV/1000000),outlier.shape = NA) +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)", color = "Fishery") +
  labs(y="",x = "Fishing PV ($Mil)") +
  coord_cartesian(xlim = c(0,500)) +
  theme_minimal()

fishsum<-dft %>% pivot_longer(cols = !`Wind farm grid ID`,names_to = "Fishery",values_to = "PV") %>% 
  group_by(Fishery) %>% 
  summarise(fishsum = sum(PV))

fishsum$Fishery<-gsub("_USD$", "",fishsum$Fishery)
fishsum$Fishery<-gsub("_", " ",fishsum$Fishery)
fishsum$Fishery<-gsub("\\.", "-",fishsum$Fishery)

bigfish<-merge(bigfish,fishsum,by="Fishery")
bigfish$PVperc<-bigfish$PV/bigfish$fishsum

bigfish45perc<-ggplot() + 
  geom_boxplot(data = bigfish, aes(y = Fishery, x = PVperc*100),outlier.shape = NA) +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)", color = "Fishery") +
  labs(y="",x = "Pecentage of total Fishing PV (%)") +
  coord_cartesian(xlim = c(0,10)) +
  theme_minimal()

rm(bigfish,dft,agree,fishroll,idim_par,selected_wfgridID,selected_rows,i,j,m,pareto_indices,s_agree,sorted_fish,fishsum)

# Plotting
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

region45<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none")

# CA 2045 ---
df<-df2[df2$`Wind farm state`=="CA",]
D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-28 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=200,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

# Plotting
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

CA45<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none")

# OR 2045 ---
df<-df2[df2$`Wind farm state`=="OR",]
D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-17 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=200,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

# Plotting
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

OR45<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none")

# WA 2045 ---
df<-df2[df2$`Wind farm state`=="WA",]
D<-as.numeric(nrow(df)) # Number of values in the function must equal D
n<-17 # Number of cells

system.time(G<-nsga2(fn=eval,
                     idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn
                     odim=2, # Output dimensions
                     lower.bounds=rep(0,nrow(df)),
                     upper.bounds=rep(1,nrow(df)),
                     popsize=200,generations=10000, cprob = 0.7, cdist = 5,
                     mprob = 0.2, mdist = 10))

# Plotting
pareto_solutions <- cbind(as.data.frame(G$value),G$pareto.optimal)
names(pareto_solutions) <- c("Sum fishing PV ($Mil)", "Mean LCOE ($/MWh)", "Pareto") # Rename columns
pareto_solutions<-pareto_solutions[pareto_solutions$Pareto=="TRUE",]
pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$`Sum fishing PV ($Mil)`), ]

WA45<-ggplot() +
  #geom_point(data = pareto_solutions, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Red"), size = 2) +
  geom_line(data = pareto_solutions_ordered, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`), color = "blue", linewidth = 1) +
  #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
  #geom_point(data = llcoe_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Black"), size = 2) +
  #geom_point(data = lfpv_pt, aes(x = `Sum fishing PV ($Mil)`, y = `Mean LCOE ($/MWh)`, color = "Green"), size = 2) +
  theme_minimal() +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
  theme(legend.title = element_blank(),legend.position="none")

# Grouping PPF and agreement plots
#region + CA + OR + WA

region45 + CA45 + OR45 + WA45 + plot_annotation(tag_levels = 'A')

(region + region45) / (CA + OR + WA) + plot_annotation(tag_levels = 'A')

bigfish30 / bigfish45 + plot_annotation(tag_levels = 'A')

bigfish30perc / bigfish45perc + plot_annotation(tag_levels = 'A')



## Scatterplot and table of LCOE per fishery
df<-read_csv("OWEP output & fishing PV data V4 DO NOT DISTRIBUTE.csv")
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

ggplot() + 
  #geom_point(data = df,aes(x = PV, y = LCOE_MWh, color = fishery)) + # Scatterplot of fisheries PV and LCOE is too jumbled with extreme values
  geom_boxplot(data = df, aes(y = fishery, x = LCOE_MWh),outlier.shape = NA) +
  #labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)", color = "Fishery") +
  labs(y="",x = "LCOE ($/MWh)") +
  theme_minimal()

## Wind results
df<-read_csv("OWEP LCOE NPV & fishing PV data V2 DO NOT DISTRIBUTE.csv")
df<-df[df$LCOE_kWh>0,] %>% drop_na(LCOE_kWh)
df$LCOE_MWh<-df$LCOE_kWh*1000
df$LCOE_kWh<-NULL
mean(df$LCOE_MWh)
df %>% group_by(State) %>% 
  summarise(mean = mean(LCOE_MWh), std = sd(LCOE_MWh))

## Region replicated across many targets
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
  #pareto_indices <-unique(as.data.frame(cbind(G$value,as.character(G$pareto.optimal))))
  pareto_indices <- which(G$pareto.optimal == TRUE)
  
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
  agree<-as.data.frame(table(as.matrix(selected_wfgridID)))
  agree$Var1<-as.numeric(levels(agree$Var1))[agree$Var1] # https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
  agree$prob<-agree$Freq/sum(agree$Freq)
  
  # Fisheries exposure at pareto selected sites (sampling approach, sample n from agreement matrix weighted by how often they are selected and calculate an expected value)
  dft<-read_csv("OWEP output & fishing PV data V4 DO NOT DISTRIBUTE.csv")
  dft<-dft %>% dplyr::select(`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD,`Weighted mean LCOE`)
  dft<-as.data.frame(dft)
  
  m<-10000 # Highest value for expected values of n
  bigfish<-data.frame("Dungeness_USD"=rep(NA,m),"At-sea_hake_USD"=rep(NA,m),"Shore_hake_USD"=rep(NA,m),"Market_squid_USD"=rep(NA,m),"Pink_shrimp_USD"=rep(NA,m),"Albacore_USD"=rep(NA,m),"Chinook_USD"=rep(NA,m),"Sablefish_USD"=rep(NA,m),"Spiny_lobster_USD"=rep(NA,m),"LCOEMWh"=rep(NA,m),"Sites"=rep(NA,m))
  for(j in 1:m){ # Sums species PV for a drawn sample of size n and iterates across samples while tracking draws in "Sites" column
    s_agree<-sample(agree$Var1,n,replace = FALSE,prob = agree$prob)
    fishroll<-data.frame("Dungeness_USD"=rep(NA,n),"At-sea_hake_USD"=rep(NA,n),"Shore_hake_USD"=rep(NA,n),"Market_squid_USD"=rep(NA,n),"Pink_shrimp_USD"=rep(NA,n),"Albacore_USD"=rep(NA,n),"Chinook_USD"=rep(NA,n),"Sablefish_USD"=rep(NA,n),"Spiny_lobster_USD"=rep(NA,n),"Weighted mean LCOE"=rep(NA,n))
    for(i in 1:n){ # Identifies species PV from a given weighted sample of size n drawn from agreement vector
      fishroll[i,]<-dft[dft$`Wind farm grid ID`==s_agree[i],2:11]
    }
    LCOEMWh<-as.data.frame(mean(fishroll$Weighted.mean.LCOE)*1000)
    fishroll$Weighted.mean.LCOE<-NULL
    fishroll<-as.data.frame(t(colSums(fishroll))) #paste(s_agree, collapse = " ")
    fishroll<-cbind(fishroll,LCOEMWh,paste(s_agree, collapse = " "))
    names(fishroll)[names(fishroll) == 'paste(s_agree, collapse = " ")']<-"Sites"
    names(fishroll)[names(fishroll) == 'mean(fishroll$Weighted.mean.LCOE) * 1000']<-"LCOEMWh"
    bigfish[j,]<-fishroll
  }
  
  LCOEMWh<-mean(bigfish$LCOEMWh)
  bigfish$LCOEMWh<-NULL
  bigfish<-bigfish %>% pivot_longer(cols = !Sites,names_to = "Fishery",values_to = "PV")
  bigfish$Fishery<-gsub("_USD$", "",bigfish$Fishery)
  bigfish$Fishery<-gsub("_", " ",bigfish$Fishery)
  bigfish$Fishery<-gsub("\\.", "-",bigfish$Fishery)
  
  bigfish<-bigfish %>% group_by(Fishery) %>% 
    summarise(meanfishPV = mean(PV))
  
  bigfish$n<-n
  bigfish$LCOEMWh<-LCOEMWh
  bigfish$dupes<-dupes
  
  return(bigfish)
}

ns<-c(seq(1,65,1))
pop<-rep(200,length(ns))

system.time(biggerfish<-map2_dfr(ns,pop,replifish))

palette <- brewer.pal(n = length(unique(biggerfish$Fishery)), name = "Paired")

ggplot(data = biggerfish, aes(x=n*.9, y=meanfishPV/1000000, group=Fishery)) + # Nominal expected impact across targets
  geom_line(aes(colour = Fishery), size = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", size = 1) +
  geom_label(aes(x = 11, y = min(meanfishPV)/1000000 - diff(range(meanfishPV)) * 0.05/1000000, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", size = 1) +
  geom_label(aes(x = 55, y = min(meanfishPV)/1000000 - diff(range(meanfishPV)) * 0.05/1000000, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Mean present value ($Mil)", colour = "Fishery")

ggplot(data = biggerfish, aes(x=n*.9, y=LCOEMWh)) + # LCOE expected impact across targets, weird artifacts from lots of pareto sites with high LCOE (TRY MEDIAN)
  geom_line(size = 1)

dft<-read_csv("OWEP output & fishing PV data V4 DO NOT DISTRIBUTE.csv")
dft<-dft %>% dplyr::select(`Wind farm grid ID`,Dungeness_USD,`At-sea_hake_USD`,Shore_hake_USD,Market_squid_USD,Pink_shrimp_USD,Albacore_USD,Chinook_USD,Sablefish_USD,Spiny_lobster_USD,`Weighted mean LCOE`)
dft<-as.data.frame(dft)
fishsum<-dft %>% pivot_longer(cols = !`Wind farm grid ID`,names_to = "Fishery",values_to = "PV") %>% 
  group_by(Fishery) %>% 
  summarise(fishsum = sum(PV))

fishsum$Fishery<-gsub("_USD$", "",fishsum$Fishery)
fishsum$Fishery<-gsub("_", " ",fishsum$Fishery)
fishsum$Fishery<-gsub("\\.", "-",fishsum$Fishery)

biggerfish<-merge(biggerfish,fishsum,by="Fishery")
biggerfish$PVperc<-biggerfish$meanfishPV/biggerfish$fishsum

ggplot(data = biggerfish, aes(x=n*.9, y=PVperc*100, group=Fishery)) + # Percent expected impact across targets
  geom_line(aes(colour = Fishery), size = 1) +
  #geom_point() + 
  geom_vline(xintercept = 11, linetype = "dashed", color = "grey50", size = 1) +
  geom_label(aes(x = 11, y = min(PVperc*100) - diff(range(PVperc*100)) * 0.05, label = "2030 Target"), fill = "white", color = "black") + # Label below the x axis
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey50", size = 1) +
  geom_label(aes(x = 55, y = min(PVperc*100) - diff(range(PVperc*100)) * 0.05, label = "2045 Target"), fill = "white", color = "black") + 
  scale_colour_manual(values = palette) +
  theme_minimal() + 
  labs(x = "Target (GW)", y = "Percent", colour = "Fishery")

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

## Sampling based approach
pareto<-function(x,n,it){
  # x - Dataframe
  # n - Subset size (# of cells needed to meet wind energy goal)
  # it - Number of sampling iterations
  
  # Storage for all evaluated solutions
  all_solutions <- data.frame(sum_v1 = numeric(), avg_v2 = numeric(), CA = numeric(), WA = numeric(), OR = numeric(), sum_MWhr_yr = numeric())
  
  for (i in 1:it) {
    # Randomly select n items
    selected_indices <- sample(nrow(x), n)
    
    # Calculate objectives
    sum_v1 <- sum(x$Total_fish_USD[selected_indices])  # Objective 1: Minimize sum of v1
    avg_v2 <- mean(x$LCOE_MWh[selected_indices])  # Objective 2: Minimize average of v2
    
    # N in each state and aggregate energy
    a<-x %>% slice(selected_indices)
    e<-sum(a$MWhyraw) # Aggregate energy
    a<-a %>% 
      group_by(State) %>% 
      summarise(count = n(),.groups = "drop") %>% 
      pivot_wider(names_from = State,values_from = count)

    # Store each solution's objectives
    all_solutions <- rbind(all_solutions, list(sum_v1 = sum_v1, avg_v2 = avg_v2, CA = a$CA, WA = a$WA, OR = a$OR, sum_MWhr_yr = e))
  }
  
  # Determine Pareto-optimality for each solution
  is_pareto_optimal <- function(index, solutions) {
    for (j in 1:nrow(solutions)) {
      if (j != index && solutions[j, "sum_v1"] <= solutions[index, "sum_v1"] && solutions[j, "avg_v2"] <= solutions[index, "avg_v2"]) {
        return(FALSE)  # Solution is dominated if another solution is better (lower) in both objectives
      }
    }
    return(TRUE)  # Solution is Pareto-optimal if no other solution is better (lower) in both objectives
  }
  
  pareto_optimal_indices <- sapply(1:nrow(all_solutions), function(i) is_pareto_optimal(i, all_solutions))
  
  # Extract Pareto-optimal solutions
  pareto_solutions <- all_solutions[pareto_optimal_indices, ]
  pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$sum_v1), ]
  
  # Adding an identifier for plotting
  all_solutions$Type <- 'Dominated'
  pareto_solutions_ordered$Type <- 'Pareto-optimal'
  
  # Combine all solutions for plotting
  plot_data <- rbind(
    all_solutions,
    pareto_solutions_ordered)
  
  a<-ggplot() +
    geom_point(data = plot_data %>% filter(Type=="Pareto-optimal"), aes(x = sum_v1, y = avg_v2, color = "Red"), size = 2) +
    geom_line(data = pareto_solutions_ordered, aes(x = sum_v1, y = avg_v2), color = "blue", linewidth = 1) +
    #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
    theme_minimal() +
    labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
    theme(legend.title = element_blank(),legend.position="none") 
  
  return(list(a,plot_data))
}

paretos<-function(x,n,it){ # Same as function above, but doesn't track outcomes by state
  # x - Dataframe
  # n - Subset size (# of cells needed to meet wind energy goal)
  # it - Number of sampling iterations
  
  # Storage for all evaluated solutions
  all_solutions <- data.frame(sum_v1 = numeric(), avg_v2 = numeric(), sum_MWhr_yr = numeric())
  
  for (i in 1:it) {
    # Randomly select n items
    selected_indices <- sample(nrow(x), n)
    
    # Calculate objectives
    sum_v1 <- sum(x$Total_fish_USD[selected_indices])  # Objective 1: Minimize sum of v1
    avg_v2 <- mean(x$LCOE_MWh[selected_indices])  # Objective 2: Minimize average of v2
    
    # Aggregate energy 
    a<-x %>% slice(selected_indices)
    e<-sum(a$MWhyraw) # Aggregate energy
    
    # Store each solution's objectives
    all_solutions <- rbind(all_solutions, list(sum_v1 = sum_v1, avg_v2 = avg_v2, sum_MWhr_yr = e))
  }
  
  # Determine Pareto-optimality for each solution
  is_pareto_optimal <- function(index, solutions) {
    for (j in 1:nrow(solutions)) {
      if (j != index && solutions[j, "sum_v1"] <= solutions[index, "sum_v1"] && solutions[j, "avg_v2"] <= solutions[index, "avg_v2"]) {
        return(FALSE)  # Solution is dominated if another solution is better (lower) in both objectives
      }
    }
    return(TRUE)  # Solution is Pareto-optimal if no other solution is better (lower) in both objectives
  }
  
  pareto_optimal_indices <- sapply(1:nrow(all_solutions), function(i) is_pareto_optimal(i, all_solutions))
  
  # Extract Pareto-optimal solutions
  pareto_solutions <- all_solutions[pareto_optimal_indices, ]
  pareto_solutions_ordered <- pareto_solutions[order(pareto_solutions$sum_v1), ]
  
  # Adding an identifier for plotting
  all_solutions$Type <- 'Dominated'
  pareto_solutions_ordered$Type <- 'Pareto-optimal'
  
  # Combine all solutions for plotting
  plot_data <- rbind(
    all_solutions,
    pareto_solutions_ordered)
  
  a<-ggplot() +
    geom_point(data = plot_data %>% filter(Type=="Pareto-optimal"), aes(x = sum_v1, y = avg_v2, color = "Red"), size = 2) +
    geom_line(data = pareto_solutions_ordered, aes(x = sum_v1, y = avg_v2), color = "blue", linewidth = 1) +
    #scale_color_manual(values = c("Pareto-optimal" = "red", "Dominated" = "gray")) +
    theme_minimal() +
    labs(x = "Sum fishing PV ($Mil)", y = "Mean LCOE ($/MWh)") +
    theme(legend.title = element_blank(),legend.position="none") 
  
  return(list(a,plot_data))
}

# Sampling
iterations<-1000
system.time(region<-pareto(x = df,n = 326, it = iterations)) # 9 hrs @ 1m iterations
system.time(ca<-paretos(x = df %>% filter(State == "CA"),n = 204, it = iterations)) # 3.5 hrs @ 1m
system.time(or<-paretos(x = df %>% filter(State == "OR"),n = 122, it = iterations)) # 3.5 hrs @ 1m
system.time(wa<-paretos(x = df %>% filter(State == "WA"),n = 122, it = iterations)) # 3.5 hrs @ 1m

(region[[1]] + ca[[1]]) / (or[[1]] + wa[[1]]) + # Unified plot
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 12))

region[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal") # Number of pareto optimal sites in each state across scenarios

region[[2]]$lcoe_x_MWh<-region[[2]]$avg_v2*region[[2]]$sum_MWhr_yr # Cost of production 
ca[[2]]$lcoe_x_MWh<-ca[[2]]$avg_v2*ca[[2]]$sum_MWhr_yr
wa[[2]]$lcoe_x_MWh<-wa[[2]]$avg_v2*wa[[2]]$sum_MWhr_yr
or[[2]]$lcoe_x_MWh<-or[[2]]$avg_v2*or[[2]]$sum_MWhr_yr

region[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal") # Number of pareto optimal sites in each state across scenarios

ca[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal") # Table of pareto optimal outcomes
or[[2]] %>% as.data.frame() %>% filter(Type == "Pareto-optimal")

df %>% group_by(State) %>% # Mean LCOE by state
  summarise(meanLCOE = mean(LCOE_MWh))

ggplot(df, aes(LCOE_MWh, colour = State)) + # CDF of LCOE by state
  stat_ecdf(linewidth = 1) +
  labs(y = "", x = "Mean LCOE ($/MWh)") +
  theme_minimal()

# ggplot(df, aes(LCOE_MWh, weight = weight, color = State)) + # Histogram by state to account for differing sample sizes
#   #geom_density() +
#   stat_bin(geom = "line", bins = 500, position = position_identity(),linewidth = 1) +
#   scale_fill_brewer(palette = "Pastel1") +
#   theme_minimal()


## GPareto - probably need to adapt the function more for our context

f<-function(x){
  print(nrow(x)) # Debug line to print the number of rows in x
  if (!is.data.frame(x)) {
    stop("The input x is not a dataframe.")
  }
  if (nrow(x) < 3) {
    stop("The dataframe x must have at least 3 rows.")
  }
  s<-x[sample(nrow(x), replace = FALSE, size = 3), ]
  sx<-sum(s$Total_fish_USD)
  sy<-mean(s$LCOE_kWh)
  
  print(s)
  print(sx)
  print(sy)
  
  return(cbind(sx,sy))
}

f2<-function(x){
  s<-as.matrix(x)
  s<-s[sample(nrow(s), replace = FALSE, size = 3), , drop = FALSE] ## See this https://stackoverflow.com/questions/50843515/sampling-columns-in-a-matrix-within-r about dropping
  sx<-sum(s[,1])
  sy<-mean(s[,2], na.rm = TRUE)
  return(cbind(sx,sy))
}

f3<-function(x){
  s<-df3[sample(nrow(df3), replace = FALSE, size = 3),]
  sx<-sum(s[,1])
  sy<-mean(s[,2], na.rm = TRUE)
  return(cbind(sx,sy))
}

f2(df3)

df3<-as.matrix(df3)

res<-easyGParetoptim(fn = f2, budget = 5, lower = c(0,0), upper = c(1,1))


####
P1 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    b1<-15*x[,1]-5
    b2<-15*x[,2]
    return(cbind((b2-5.1*(b1/(2*pi))^2+5/pi*b1-6)^2 +10*((1-1/(8*pi))*cos(b1)+1),
                 -sqrt((10.5-b1)*(b1+5.5)*(b2+0.5)) - 1/30*(b2 -5.1*(b1/(2*pi))^2-6)^2 - 1/3*((1-1/(8*pi))*cos(b1)+1)
    ) 
    )
  }

d <- 2; ninit <- 10; fun <- P1
design <- DiceDesign::lhsDesign(ninit, d, seed = 42)$design
class(design)

res<-easyGParetoptim(fn = P1, budget = 50, lower = c(0,0), upper = c(1,1))
plotGPareto(res)

#####
MOP2 <- function(x)
{
  xmod <- x*4 - 2
  if (is.null(nrow(x)))
  { 
    n <- length(xmod)
    y1 <- 1 - exp(-sum((xmod - 1/sqrt(n))^2) )
    y2 <- 1 - exp(-sum((xmod + 1/sqrt(n))^2) )
    Y <- matrix(c(y1,y2),1,2)
  } else
  {
    n <- ncol(xmod)
    y1 <- 1 - exp(-rowSums((xmod - 1/sqrt(n))^2) )
    y2 <- 1 - exp(-rowSums((xmod + 1/sqrt(n))^2) )
    Y <- cbind(y1,y2)
  }
  
  return(Y)
}

xmod <- a*4 - 2

design.init <- matrix(seq(0, 1, length.out = 6), ncol = 1)
response.init <- MOP2(design.init)
MOP2(design.init)

mf1 <- km(~1, design = design.init, response = response.init[, 1])
mf2 <- km(~1, design = design.init, response = response.init[, 2])
model <- list(mf1, mf2)

res <- GParetoptim(model = model, fn = MOP2, crit = "EHI", nsteps = 25,lower = 0, upper = 1, critcontrol = list(refPoint = c(2, 2)))
plotGPareto(res)
res<-easyGParetoptim(fn = MOP2, budget = 100, lower = c(0,0), upper = c(1,1))
plotGPareto(res)


###############

DTLZ2 <- function(x, nobj = 3){
  if(is.null(dim(x))){
    x <- matrix(x, 1) 
  }
  n <- ncol(x)
  
  y <- matrix(x[,1:(nobj-1)], nrow(x))
  z <- matrix(x[,nobj:n], nrow(x))
  
  g <- rowSums((z-0.5)^2)
  
  #   tmp <- c(rev(cumprod(cos(y * pi/2))), 1)
  #   tmp2 <- c(1, rev(sin(y * pi/2)))
  tmp <- t(apply(cos(y * pi/2), 1, cumprod))
  tmp <- cbind(t(apply(tmp, 1, rev)), 1)
  
  tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
  
  f <- tmp * tmp2 * (1 + g)
  
}

y <- matrix(data = x[,1:(3-1)], nrow(x))
z <- matrix(x[,3:n], nrow(x))

a<-rbind(rep(1,5),seq(1,5,1))



DTLZ2

ZDT2(a)

?ZDT3()

    
