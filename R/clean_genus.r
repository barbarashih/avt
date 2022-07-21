#' Remove poorly defined Genus from the taxonomy table
#'
#' Clean up Genus to remove poorly defined genus
#' 
#' @param x Taxonomy dataframe.
#' @param y A vector with of Genus that should be ignored.

#' @export
cleanGenus <- function(x, y = vector()){
	default_genus_to_ignore = c("Allorhizobium.Neorhizobium", 
								"bacterium episymbiont of Kiwa sp.",
								"Bacteroidia bacterium canine oral taxon 187",
								"CAG-352",
								"CL500-3",
								"Clade Ia",
								"F0058",
								"F0332",
								"Family XIII UCG-001",
								"IS-44",
								"Labrys",
								"Mastigocladopsis PCC-10914",
								"MB11C04 marine group",
								"NS10 marine group",
								"NS2b marine group",
								"NS3a marine group",
								"NS4 marine group",
								"NS5 marine group",
								"OM43 clade",
								"OM60(NOR5) clade",
								"Pir4 lineage" ,
								"PMMR1",
								"R76-B128",
								"RB41" ,
								"RS62 marine group",
								"SAR92 clade",
								"Schizothrix LEGE 07164",
								"SM1A02",
								"uncultured bacterium",
								"uncultured Cytophagales bacterium",
								"uncultured marine bacterium",
								"uncultured organism",
								"uncultured PS1 clade bacterium",
								"uncultured soil bacterium",
								"uncultured Verrucomicrobia bacterium",
								"uncultured",
								"unidentified",
								"Verruc-01",
								"W5053",
								"WCHB1-32",
								"hydrothermal vent metagenome",
								"Acrophormium PCC-7375",
								"CHKCI002",
								"environmental clone OCS162",
								"environmental clone OCS182",
								"Sva0996 marine group",
								"UC5-1-2E3",
								"Ulva sp. UNA00071828",
								"uncultured bacterium adhufec202",
								"uncultured bacterium ARCTIC10_F_05",
								"uncultured Lentisphaerae bacterium",
								"uncultured marine microorganism",
								"uncultured methanogenic archaeon",
								"uncultured Thermoanaerobacterales bacterium",
								"unidentified rumen bacterium RF39",
								"unidentified marine bacterioplankton",
								"CHKCI001")
	for(ignore_genus in default_genus_to_ignore){
		x<-x[!x$Genus %in% ignore_genus,]
	}

	
	x$Genus <- gsub("^uncultured Bacteroidetes.*?$", "Bacteroidetes bacterium", x$Genus)
	x$Genus <- gsub("^uncultured Bacteroidales.*?$", "Bacteroidales bacterium", x$Genus)
	x$Genus <- gsub("^uncultured Flavobacteriales.*?$", "Flavobacteria bacterium", x$Genus)
	x$Genus <- gsub("^uncultured alpha proteobacterium.*?$", "Proteobacterium bacterium", x$Genus)
	x$Genus <- gsub("^uncultured Sphingomonas sp.*?$", "Sphingomonas", x$Genus)
	x$Genus <- gsub("^uncultured Thermoanaerobacterales*?$", "Thermoanaerobacterales bacterium", x$Genus)
	x$Genus <- gsub("^uncultured Mollicutes bacterium*?$", "Mollicutes bacterium", x$Genus)
	x$Genus <- gsub("^uncultured eubacterium bacterium*?$", "Eubacterium", x$Genus)
	x$Genus <- gsub("^Alphaproteobacteria bacterium.*?$", "Proteobacterium bacterium", x$Genus)
	x$Genus <- gsub("^Clostridiales bacterium.*?$", "Clostridiales bacterium", x$Genus)
	x$Genus <- gsub("^Coriobacteriaceae.*?$", "Coriobacteriaceae bacterium", x$Genus) # unclassified
	x$Genus <- gsub("^Proteobacteria bacterium.*?$", "Proteobacteria bacterium", x$Genus) # unclassified
	x$Genus <- gsub("^Azospirillum sp. 47_25.*?$", "Azospirillum", x$Genus) # unclassified
	x$Genus <- gsub("^candidate division SR1 bacterium MGEHA.*?$", "Candidatus Absconditabacteria", x$Genus) # unclassified
	x$Genus <- gsub("^Gracilibacteria bacterium oral taxon 873.*?$", "Candidatus Gracilibacteria", x$Genus) # unclassified
	x$Genus <- gsub("^S5-A14a.*?$", "Clostridiales bacterium", x$Genus) # unclassified
	x$Genus <- gsub("^GCA-900066575.*?$", "Clostridiales bacterium", x$Genus) # unclassified
	x$Genus <- gsub("^TM7 phylum sp. oral clone FR058.*?$", "Candidatus Saccharibacteria bacterium", x$Genus) # unclassified
	x$Genus <- gsub("^DTU089.*?$", "Clostridiales bacterium", x$Genus) # unclassified
	x$Genus <- gsub("^Family XIII AD3011 group.*?$", "Clostridiales bacterium", x$Genus) # unclassified
	x$Genus <- gsub("^uncultured eubacterium.*?$", "Eubacterium", x$Genus) # unclassified
	x$Genus <- gsub("^gut metagenome*?$", "Other", x$Genus)
	x$Genus <- gsub("^Prevotellaceae .*?$", "Prevotella", x$Genus)
	x$Genus <- gsub("^Ruminococcaceae .*?$", "Ruminococcus", x$Genus)
	x$Genus <- gsub("^Erysipelotrichaceae .*?$", "Erysipelotrichia", x$Genus)
	x$Genus <- gsub("^Lachnospiraceae .*?$", "Lachnospira", x$Genus)
	x$Genus <- gsub("^Synechococcus .*?$", "Synechococcus", x$Genus)
	x$Genus <- gsub("^Rikenellaceae .*?$", "Rikenella", x$Genus)
	x$Genus <- gsub("^Rivularia .*?$", "Rivularia", x$Genus)
	x$Genus <- gsub("^Sphingomonas .*?$", "Sphingomonas", x$Genus)
	x$Genus <- gsub("^Christensenellaceae.*?$", "Christensenella", x$Genus)
	x$Genus <- gsub("^Porphyromonadaceae.*?$", "Porphyromonas", x$Genus)
	x$Genus <- gsub("^Pleurocapsa.*?$", "Pleurocapsa", x$Genus)
	x$Genus <- gsub("^Dehalococcoidia.*?$", "Dehalococcus", x$Genus)
	x$Genus <- gsub("^Clostridium.*?$", "Clostridium", x$Genus)
	x$Genus <- gsub("^Defluviitaleaceae.*?$", "Defluviitaleaceae", x$Genus)
	x$Genus <- gsub("^Atelocyanobacterium.*?$", "Atelocyanobacterium", x$Genus)
	x$Genus <- gsub("^Musa acuminata.*?$", "Musa acuminata", x$Genus)
	x$Genus <- gsub("^Pleurocapsa.*?$", "Pleurocapsa", x$Genus)
	x$Genus <- gsub("^Wilmottia.*?$", "Wilmottia", x$Genus)
	x$Genus <- gsub("^GCA-900066755$", "Clostridium", x$Genus)
	x$Genus <- gsub("^uncultured Clostridia bacterium$", "Clostridia bacterium", x$Genus)
	x$Genus <- gsub("^uncultured Clostridia bacterium$", "Clostridia bacterium", x$Genus)
	x$Genus <- gsub("^uncultured crenarchaeote$", "crenarchaeote", x$Genus)
	x$Genus <- gsub("^uncultured cyanobacterium$", "cyanobacterium", x$Genus)
	x$Genus <- gsub("\\].*?$", "", x$Genus)
	x$Genus <- gsub("^\\[", "", x$Genus)
	x$Genus <- gsub("^Gracilibacteria bacterium.*?$", "Gracilibacteria bacterium", x$Genus)


	return(x)
}


