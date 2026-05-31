.PHONY: help test test-fast lint dryrun report clean dag fetch-test-data

SNAKEMAKE ?= snakemake
CORES     ?= 4
PROFILE   ?=

help:
	@echo "germline-flow targets:"
	@echo "  make test          run the bundled chr20 end-to-end test (uses test/config.yaml)"
	@echo "  make test-fast     same as 'test' but reuses any existing intermediate files"
	@echo "  make lint          run snakemake --lint on the workflow"
	@echo "  make dryrun        snakemake -n on the default config"
	@echo "  make dag           render the DAG as docs/dag.svg (requires graphviz)"
	@echo "  make report        rebuild only the HTML report rule"
	@echo "  make clean         remove results/ and .snakemake/ (does NOT touch resources/)"
	@echo ""
	@echo "Variables:"
	@echo "  CORES=$(CORES)  SNAKEMAKE=$(SNAKEMAKE)  PROFILE=$(PROFILE)"

# ---- linting -----------------------------------------------------------------

lint:
	$(SNAKEMAKE) --lint --configfile config/config.yaml

# ---- dry run -----------------------------------------------------------------

dryrun:
	$(SNAKEMAKE) -n --configfile config/config.yaml --use-conda

# ---- end-to-end test ---------------------------------------------------------
# The test config uses a downsampled chr20 reference and three tiny FASTQ pairs
# bundled under test/. CI runs exactly this target.

test:
	$(SNAKEMAKE) \
	    --configfile test/config.yaml \
	    --use-conda --conda-frontend mamba \
	    --cores $(CORES) \
	    --rerun-incomplete \
	    --show-failed-logs \
	    all

test-fast:
	$(SNAKEMAKE) \
	    --configfile test/config.yaml \
	    --use-conda --conda-frontend mamba \
	    --cores $(CORES) \
	    --keep-going \
	    all

# ---- report only -------------------------------------------------------------

report:
	$(SNAKEMAKE) \
	    --configfile config/config.yaml \
	    --use-conda \
	    --cores 1 \
	    --forcerun report_html \
	    results/report/cohort_qc.html

# ---- DAG visualisation -------------------------------------------------------

dag:
	mkdir -p docs
	$(SNAKEMAKE) --configfile test/config.yaml --dag all | dot -Tsvg > docs/dag.svg
	@echo "wrote docs/dag.svg"

# ---- cleanup -----------------------------------------------------------------

clean:
	rm -rf results/ .snakemake/
	@echo "removed results/ and .snakemake/ (resources/ left intact)"
