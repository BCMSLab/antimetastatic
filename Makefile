#!/bin/bash
DATA=data
LOG=log
SCRIPTS=scripts

all: directories \
	$(LOG)/clean_data.R.Rout \
	$(LOG)/knocktf.R.Rout \
	$(LOG)/hms.R.Rout \
	$(LOG)/lincs.R.Rout \
	$(LOG)/lincs_diff_expr.R.Rout \
	$(LOG)/diff_expr_groups.R.Rout \
	$(LOG)/npa_model.R.Rout \
	$(LOG)/scoring.R.Rout \
	$(LOG)/concordance.R.Rout \
	$(LOG)/coherence.R.Rout \
	manuscript
	
directories:
	mkdir -p $(DATA)
	mkdir -p $(LOG)
	@echo "Make two directories data/ and log/."
$(LOG)/clean_data.R.Rout: $(SCRIPTS)/clean_data.R \
	$(DATA)/md.csv
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Download and clean gene expression data."
$(LOG)/diff_expr_groups.R.Rout: $(SCRIPTS)/diff_expr_groups.R \
	$(DATA)/msg_kd_eset.rds \
	$(DATA)/msg_list.txt
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Perform differential expression analysis."
$(LOG)/knocktf.R.Rout: $(SCRIPTS)/knocktf.R
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Download KnockTF data."
$(LOG)/hms.R.Rout: $(SCRIPTS)/hms.R \
	$(DATA)/dataset_20268_20210412113444.csv \
	$(DATA)/dataset_20269_20210412114601.csv
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Download HMS data."	
$(LOG)/lincs.R.Rout: $(SCRIPTS)/lincs.R
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Download LINCS data."
$(LOG)/lincs_diff_expr.R.Rout: $(SCRIPTS)/lincs_diff_expr.R \
	$(DATA)/lincs_trt_cp.rds
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Perfom differential expression on LINCS data."
$(LOG)/npa_model.R.Rout: $(SCRIPTS)/npa_model.R \
	$(DATA)/combined_interactions.csv \
	$(DATA)/diff_expr_groups.rds \
	$(DATA)/knock_tf_subset.tsv
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Format the network for NPAModels"
$(LOG)/scoring.R.Rout: $(SCRIPTS)/scoring.R \
	$(DATA)/lincs_trt_cp_deg.rds \
	$(DATA)/full_graph.rds \
	$(DATA)/lincs_trt_cp_deg.tsv \
	$(LOG)/npa_model.R.Rout
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Score the network perturbations."
$(LOG)/concordance.R.Rout: $(SCRIPTS)/concordance.R \
	$(DATA)/combined_interactions.csv \
	$(LOG)/npa_model.R.Rout \
	$(LOG)/scoring.R.Rout
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Calculate network concordance."
$(LOG)/coherence.R.Rout: $(SCRIPTS)/coherence.R \
	$(DATA)/combined_interactions.csv \
	$(LOG)/npa_model.R.Rout \
	$(LOG)/scoring.R.Rout
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Calculate network coherence."
manuscript: manuscript.Rmd \
	$(LOG)/npa_model.R.Rout \
	$(LOG)/scoring.R.Rout
	Rscript -e "knitr::knit('manuscript.Rmd')"
	