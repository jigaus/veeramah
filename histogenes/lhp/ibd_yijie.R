library(dplyr)
library(stringr)
library(igraph)
library(tidyr)

mcs_ex <- read.csv("mcs_ext.csv", header = T)

### FILTERING PEDIGREES -------------------------------------------------------------------------

# sites
sites <- c("rko", "rkf", "hnj", "tap", "led", "rkc", "mgs", "kup", "kfp", "csk", "mdf", "leo", "tgh", "kih", "ptl", "nja")
sites <- sort(sites)

# filtering
rko <- dplyr::filter(mcs_ex, grepl("RKO", mcs_ex$iid2))
rkf <- dplyr::filter(mcs_ex, grepl("RKF", mcs_ex$iid2))
hnj <- dplyr::filter(mcs_ex, grepl("HNJ", mcs_ex$iid2))
tap <- dplyr::filter(mcs_ex, grepl("Tap", mcs_ex$iid2))
led <- dplyr::filter(mcs_ex, grepl("Ledine", mcs_ex$iid2))
rkc <- dplyr::filter(mcs_ex, grepl("RKC", mcs_ex$iid2))
mgs <- dplyr::filter(mcs_ex, grepl("MGS", mcs_ex$iid1))
kup <- dplyr::filter(mcs_ex, grepl("KUP", mcs_ex$iid1))
kfp <- dplyr::filter(mcs_ex, grepl("KFP", mcs_ex$iid1))
csk <- dplyr::filter(mcs_ex, grepl("CSK", mcs_ex$iid1))
mdf <- dplyr::filter(mcs_ex, grepl("MDF", mcs_ex$iid1))
leo <- dplyr::filter(mcs_ex, grepl("LEO", mcs_ex$iid1))
tgh <- dplyr::filter(mcs_ex, grepl("TGH", mcs_ex$iid1))
kih <- dplyr::filter(mcs_ex, grepl("KIH", mcs_ex$iid1))
ptl <- dplyr::filter(mcs_ex, grepl("PTL", mcs_ex$iid1))
nja <- dplyr::filter(mcs_ex, grepl("NJA", mcs_ex$iid1))

# connections
connections <- data.frame(
  sites = sites,
  number = c(71, 3, 16, 33, 7, 2, 14, 36, 459, 4, 14, 19, 22, 20, 1, 18)
)
sum(connections$number) # check lol

barplot(
  height = connections$number,
  names.arg = connections$sites,
  xlab = "sites",
  ylab = "number of ibd connections"
)

## pedigree iii
three <- c(644, 748, 647, 648, 671, 645, 744, 701, 708, 736, 705, 711, 759, 751)
three_connect_one <- mcs_ex %>% filter(as.numeric(str_extract(iid1, "\\d+")) %in% three)
three_connect_two <- mcs_ex %>% filter(as.numeric(str_extract(iid2, "\\d+")) %in% three)
three_connect <- rbind(three_connect_one, three_connect_two)

iii_connect <- data.frame(
  source = three_connect$iid1,
  target = three_connect$iid2
)
iii_network <- graph_from_data_frame(iii_connect, directed = F)
plot(iii_network)

write.table(three_connect, append = F, file = "three_connect.csv", quote = F, sep = ",", row.names = F)

## pedigree ii 
two <- c(553, 539, 543, 643, 566, 688, 596, 682, 642, 694, 696, 697, 658, 556, 693, 571, 659, 571, 659, 552, 517)
two_connect_one <- mcs_ex %>% filter(as.numeric(str_extract(iid1, "\\d+")) %in% two)
two_connect_two <- mcs_ex %>% filter(as.numeric(str_extract(iid2, "\\d+")) %in% two)
two_connect <- rbind(two_connect_one, two_connect_two)

write.table(two_connect, file = "two_connect.csv", quote = F, sep = ",", row.names = F)

ii_connect <- data.frame(
  source = two_connect$iid1,
  target = two_connect$iid2
)
ii_network <- graph_from_data_frame(ii_connect, directed = F)
plot(ii_network,
     vertex.size = 20)


plot(ii_network, vertex.size = 10, 
vertex.label.color = "black", vertex.label.cex = 0.8, vertex.label.degree = -pi/2,
edge.arrow.size = 0.7, edge.arrow.width = 0.4, edge.color = "black") 

## pedigree i 
one <- c(700, 704, 684, 432, 454, 703, 702, 491, 450, 732, 637, 430, 435, 455, 484, 433, 533, 534, 458, 489, 490, 563, 456, 600, 599, 486, 602, 459, 449)
one_connect_one <- mcs_ex %>% filter(as.numeric(str_extract(iid1, "\\d+")) %in% one)
one_connect_two <- mcs_ex %>% filter(as.numeric(str_extract(iid2, "\\d+")) %in% one)
one_connect <- rbind(one_connect_one, one_connect_two)

i_connect <- data.frame(
  source = one_connect$iid1,
  target = one_connect$iid2
)
i_network <- graph_from_data_frame(i_connect, directed = F)
plot(i_network)

plot(ii_network, vertex.size = 9, vertex.color = rainbow(10, .8, .8, alpha = .8),
vertex.label.color = "black", vertex.label.cex = 0.8, vertex.label.degree = -pi/4,
edge.arrow.size = 0.4, edge.arrow.width = 0.6, edge.color = "black") 

write.table(one_connect, append = F, file = "one_connect.csv", quote = F, sep = ",", row.names = F)

## small pedigrees
small <- c(630, 471, 473, 466, 472, 620, 635, 716, 714, 715, 522, 520, 545, 570, 540, 542, 535, 575, 650, 453, 464)
small_connect_one <- mcs_ex %>% filter(as.numeric(str_extract(iid1, "\\d+")) %in% small)
small_connect_two <- mcs_ex %>% filter(as.numeric(str_extract(iid2, "\\d+")) %in% small)
small_connect <- rbind(small_connect_one, small_connect_two)

s_connect <- data.frame(
  source = small_connect$iid1,
  target = small_connect$iid2
)
s_network <- graph_from_data_frame(s_connect, directed = F)
plot(
  s_network,
  edge.color = "black"
)

write.table(small_connect, append = F, file = "small_connect.csv", quote = F, sep = ",", row.names = F)

### STATISTICAL TESTS ------------------------------------------------
# as stringent as possible when considering IBD > 12
# removing all individuals not considered for sites of interest
# only considering individuals in pedigree with traceable ancestry

admix <- read.csv("7thCEsamples_all_Baltic.csv", header = T) # new, with mcs
admix <- read.csv("admix.csv", header = T) # old 
admix <- subset(admix, select = -c(Index))

mcs_list <- read.csv("mcs_list.csv", header = F)
mcs_list <- as.vector(mcs_list)

# vector of unique indiv across pedigrees
# filter admix tables by vector(s)

## pedigree iii
three_connect <- read.csv("three_connect.csv", header = T) # removed any PTL, reload above for og 
three_connect <- subset(three_connect, select = -c(n_IBD_12, sum_IBD_20, n_IBD_20)) # removing columns for looking
out_iii <- as.vector(three_connect$iid1)
in_iii <- as.vector(three_connect$iid2)

out_iii_admix <- admix[match(out_iii, admix$ID), ]
colnames(out_iii_admix)[2] <- "CASIA_out"
colnames(out_iii_admix)[3] <- "EASIA_out"
colnames(out_iii_admix)[4] <- "GRK_out"
colnames(out_iii_admix)[5] <- "MEDEU_out"
colnames(out_iii_admix)[6] <- "NAFRICA_out"
colnames(out_iii_admix)[7] <- "NGBI_out"
colnames(out_iii_admix)[8] <- "SASIA_out"
colnames(out_iii_admix)[9] <- "SCAND_out"
colnames(out_iii_admix)[10] <- "SUBSAHARAN_out"

in_iii_admix <- admix[match(in_iii, admix$ID), ]
colnames(in_iii_admix)[2] <- "CASIA_in"
colnames(in_iii_admix)[3] <- "EASIA_in"
colnames(in_iii_admix)[4] <- "GRK_in"
colnames(in_iii_admix)[5] <- "MEDEU_in"
colnames(in_iii_admix)[6] <- "NAFRICA_in"
colnames(in_iii_admix)[7] <- "NGBI_in"
colnames(in_iii_admix)[8] <- "SASIA_in"
colnames(in_iii_admix)[9] <- "SCAND_in"
colnames(in_iii_admix)[10] <- "SUBSAHARAN_in"

three_total <- cbind(three_connect, in_iii_admix)
three_total <- subset(three_total, select = -c(ID))
three_total <- cbind(three_connect, out_iii_admix)
three_total <- subset(three_total, select = -c(ID))

plot(
  x = three_total$sum_IBD_12,
  y = three_total$NGBI_out
)

## pedigree ii 
two_connect <- read.csv("two_connect.csv", header = T) # removed any 696B, reload above for og 
two_connect <- subset(two_connect, select = -c(n_IBD_12, sum_IBD_20, n_IBD_20)) # removing columns for looking
out_ii <- as.vector(two_connect$iid1)
in_ii <- as.vector(two_connect$iid2)

out_ii_admix <- admix[match(out_ii, admix$ID), ]
colnames(out_ii_admix)[2] <- "CASIA_out"
colnames(out_ii_admix)[3] <- "EASIA_out"
colnames(out_ii_admix)[4] <- "GRK_out"
colnames(out_ii_admix)[5] <- "MEDEU_out"
colnames(out_ii_admix)[6] <- "NAFRICA_out"
colnames(out_ii_admix)[7] <- "NGBI_out"
colnames(out_ii_admix)[8] <- "SASIA_out"
colnames(out_ii_admix)[9] <- "SCAND_out"
colnames(out_ii_admix)[10] <- "SUBSAHARAN_out"

in_ii_admix <- admix[match(in_ii, admix$ID), ]
colnames(in_ii_admix)[2] <- "CASIA_in"
colnames(in_ii_admix)[3] <- "EASIA_in"
colnames(in_ii_admix)[4] <- "GRK_in"
colnames(in_ii_admix)[5] <- "MEDEU_in"
colnames(in_ii_admix)[6] <- "NAFRICA_in"
colnames(in_ii_admix)[7] <- "NGBI_in"
colnames(in_ii_admix)[8] <- "SASIA_in"
colnames(in_ii_admix)[9] <- "SCAND_in"
colnames(in_ii_admix)[10] <- "SUBSAHARAN_in"

two_total <- cbind(two_connect, in_ii_admix)
two_total <- subset(two_total, select = -c(site))
two_total <- cbind(two_connect, out_ii_admix)
two_total <- subset(two_total, select = -c(ID))

## pedigree i
one_connect <- read.csv("one_connect.csv", header = T) # removed PTL, TGH, Tap, reload above for OG

## the big table 
# up -> down = 1, 2, 3

data <- data.frame(
  ibd_connections = c(191, 72, 27),
  sample_size = c(29, 19, 14)
)

plot(
  x = data$sample_size,
  y = data$ibd_connections,
  xlab = "sample size",
  ylab = "number of ibd counts/pedigree"
)

poi <- glm(
  ibd_connections ~ sample_size,
  data = data,
  family = poisson(link = "log")
)

summary(poi)

# divide abs number / # of possible connections = normalization 
# raw / possible
