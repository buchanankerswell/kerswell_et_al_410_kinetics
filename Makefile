# Top-level dirs
PROJECT_ROOT := $(CURDIR)
R := $(PROJECT_ROOT)/R
BASH := $(PROJECT_ROOT)/bash
PYTHON := $(PROJECT_ROOT)/python
SIMULATION := $(PROJECT_ROOT)/simulation
DRAFT := $(PROJECT_ROOT)/draft

# Targets
.PHONY: install uninstall models plume-models slab-models adiabat sync-barkla download visualize manuscript environments clean deep-clean help

build: environments download adiabat visualize
	@open figs/*.png

install:
	@$(MAKE) --no-print-directory -C $(BASH)

uninstall:
	@$(MAKE) --no-print-directory -C $(BASH) uninstall

models:
	@$(MAKE) --no-print-directory -C $(SIMULATION) models

plume-models:
	@$(MAKE) --no-print-directory -C $(SIMULATION) plume-models

slab-models:
	@$(MAKE) --no-print-directory -C $(SIMULATION) slab-models

adiabat:
	@$(MAKE) --no-print-directory -C $(PYTHON) adiabat

sync-barkla:
	@$(MAKE) --no-print-directory -C $(PYTHON) sync-barkla

download:
	@$(MAKE) --no-print-directory -C $(R) download

visualize:
	@$(MAKE) --no-print-directory -C $(R) visualize
	@$(MAKE) --no-print-directory -C $(PYTHON) visualize

manuscript:
	@$(MAKE) --no-print-directory -C $(DRAFT)

environments:
	@echo "    --------------------------------------------------"
	@echo "    Creating conda environments"
	@echo "    --------------------------------------------------"
	@if conda info --envs | awk '{print $$1}' | grep -qx "kerswell-et-al-410-kinetics"; then \
	  echo " -- Found conda environment: kerswell-et-al-410-kinetics"; \
	else \
	  echo " -> Creating conda environment: kerswell-et-al-410-kinetics from environment.yaml"; \
	  conda env create -f $(PYTHON)/environment.yaml; \
	fi
	@echo "    --------------------------------------------------"
	@echo "    Creating R environment"
	@echo "    --------------------------------------------------"
	@Rscript $(R)/environment.R

clean:
	@echo "    --------------------------------------------------"
	@echo "    Cleaning ..."
	@echo "    --------------------------------------------------"
	@$(MAKE) --no-print-directory -C $(BASH) clean || true
	@$(MAKE) --no-print-directory -C $(PYTHON) clean || true
	@$(MAKE) --no-print-directory -C $(R) clean || true
	@$(MAKE) --no-print-directory -C $(SIMULATION) clean || true
	@$(MAKE) --no-print-directory -C $(DRAFT) clean || true
	@find . -name ".DS_Store" -type f -delete

deep-clean: clean

	@echo "    --------------------------------------------------"
	@echo "    Deep cleaning ..."
	@echo "    --------------------------------------------------"
	@$(MAKE) --no-print-directory -C $(BASH) deep-clean || true
	@$(MAKE) --no-print-directory -C $(PYTHON) deep-clean || true
	@$(MAKE) --no-print-directory -C $(R) deep-clean || true
	@$(MAKE) --no-print-directory -C $(SIMULATION) deep-clean || true
	@$(MAKE) --no-print-directory -C $(DRAFT) deep-clean || true

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    install       Build deal.II and ASPECT"
	@echo "    uninstall     Uninstall deal.II and ASPECT"
	@echo "    models        Run ASPECT 2d box models"
	@echo "    plume-models  Run ASPECT 2d box plume models"
	@echo "    slab-models   Run ASPECT 2d box slab models"
	@echo "    adiabat       Generate adiabatic profiles for ASPECT models"
	@echo "    sync-barkla   Sync data from barkla2 cluster (UoL)"
	@echo "    download      Download ASPECT results from OSF repo"
	@echo "    visualize     Visualize all results"
	@echo "    environments  Create Conda environments"
	@echo "    clean         Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean    Deep clean figures and results (use with caution!)"
	@echo "    help          Show this help message"
