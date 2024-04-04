README - meta data



01_depth.R
	goal: create graph of sampling depths between eDNA + trawl 
	
	inputs: 
	"trawl_metadata.csv"
	"eDNA_metadata.csv"
	"lat_lon_all.csv"
	
	outputs: 
	"samplingdepths_all_.png"
	
02_maps.R
	goal: create map of study site 
	
	inputs: 
	"trawl_metadata.csv"
	"eDNA_metadata.csv"
	"lat_lon_all.csv"
	
	outputs: 
	mapbw.png
	

03_ASV_reads.R
	goal: create datafile that contains combined reads in, reads out and percentage
	retained for 12se and 12su data
	
	inputs: 
	retained_reads.filterAndTrim_step.length_var.txt (for 12se + 12su )
	
	outputs: 
	reads_retained.csv