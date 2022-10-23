library(phyloseq)
library(microbiome)
library(RColorBrewer)


otu.tab <- read_tsv("https://raw.githubusercontent.com/krabberod/UNIS_AB332_2022/main/computer_lab/data/AB332_otutab_reduc3.txt")
otu.tab <- column_to_rownames(otu.tab, var = "OTUNumber")

tax.tab<-read_tsv("https://raw.githubusercontent.com/krabberod/UNIS_AB332_2022/main/computer_lab/data/AB332_2021_taxtab.txt", col_names=F)
colnames(tax.tab) <- c("OTUname","acc","Kingdom","Supergroup","Phylum","Tax1",
                       "Tax2","Species","Id","E-val","bit","length")
tax.tab<-column_to_rownames(tax.tab,var="OTUname")
OTU<-otu_table(otu.tab,taxa_are_rows=T)
TAX<-tax.tab %>% dplyr::select(Kingdom:Species) %>% as.matrix %>% tax_table()

isa.phyloseq<-phyloseq(OTU,TAX)


# Calculate the relative abundance: 
isa.phyloseq.rel <- transform(isa.phyloseq, "compositional")

# Example: select the Dinoflagelleates
isa.dinos <- subset_taxa(isa.phyloseq.rel, Phylum %in% "Dinoflagellates")

# Example: select MALV
isa.MALV <- subset_taxa(isa.phyloseq.rel, Phylum %in% "MALVs")

# Example: select micromonas
isa.micromonas <- subset_taxa(isa.phyloseq.rel, Species %in% c("Micromonas_CCMP2099_Arctic","Micromonas_CCMP1195_clade_C","Micromonas_CCMP1545_clade_D"))

# Plotting the entire community
plot_bar(isa.phyloseq.rel,fill = "Supergroup") +
  scale_fill_manual(values = sample(colorRampPalette(brewer.pal(8,"BrBG"))(14)),name="Supergroup") +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 8)+
  theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.title.x = element_blank() ,axis.text.x  = element_text(angle=90)) +
  ylab("Relative Abundance\n") +
  ggtitle("All groups")

# Plotting the Malvs: 
plot_bar(isa.MALV,fill = "Species") +
  scale_fill_manual(values = sample(colorRampPalette(brewer.pal(8,"BrBG"))(21)),name="Species") +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 8)+
  theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.title.x = element_blank() ,axis.text.x  = element_text(angle=90,)) +
  ylab("Relative Abundance\n") +
  ggtitle("MALV")

# Plotting the Micromonas: 
plot_bar(isa.micromonas,fill = "Species") +
  scale_fill_manual(values = sample(colorRampPalette(brewer.pal(8,"BrBG"))(3)),name="Species") +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 8)+
  theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.title.x = element_blank() ,axis.text.x  = element_text(angle=90)) +
  ylab("Relative Abundance\n") +
  ggtitle("Micromonas")

