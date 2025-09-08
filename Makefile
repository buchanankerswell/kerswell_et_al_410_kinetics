# Top-level dirs
PROJECT_ROOT := $(CURDIR)
R := $(PROJECT_ROOT)/R
BASH := $(PROJECT_ROOT)/bash
PYTHON := $(PROJECT_ROOT)/python
SIMULATION := $(PROJECT_ROOT)/simulation
MARKDOWN := $(PROJECT_ROOT)/markdown

# Targets
.PHONY: install uninstall 2d-box-models 2d-shell-models adiabat visualize environments clean deep-clean help

install:
	@$(MAKE) --no-print-directory -C $(BASH)

uninstall:
	@$(MAKE) --no-print-directory -C $(BASH) uninstall

2d-shell-models:
	@$(MAKE) --no-print-directory -C $(SIMULATION) 2d-shell-models

2d-box-models:
	@$(MAKE) --no-print-directory -C $(SIMULATION) 2d-box-models

adiabat:
	@$(MAKE) --no-print-directory -C $(PYTHON) adiabat

sync-barkla:
	@$(MAKE) --no-print-directory -C $(PYTHON) sync-barkla

visualize:
	@$(MAKE) --no-print-directory -C $(R) visualize
	@$(MAKE) --no-print-directory -C $(PYTHON) visualize

manuscript-diff:
	@$(MAKE) --no-print-directory -C $(MARKDOWN) diff

manuscript:
	@$(MAKE) --no-print-directory -C $(MARKDOWN)

environments:
	@for file in $(PYTHON)/environments/*.yml; do \
		name=$$(basename $$file .yml); \
		if conda info --envs | awk '{print $$1}' | grep -qx "$$name"; then \
			echo " -- Environment '$$name' already exists. Skipping ..."; \
		else \
			echo " -> Creating environment: $$name from $$file"; \
			conda env create -n "$$name" -f "$$file"; \
		fi \
	done

clean:
	@echo "    --------------------------------------------------"
	@echo "    Cleaning ..."
	@echo "    --------------------------------------------------"
	@$(MAKE) --no-print-directory -C $(BASH) clean || true
	@$(MAKE) --no-print-directory -C $(PYTHON) clean || true
	@$(MAKE) --no-print-directory -C $(R) clean || true
	@$(MAKE) --no-print-directory -C $(SIMULATION) clean || true
	@$(MAKE) --no-print-directory -C $(MARKDOWN) clean || true
	@find . -name ".DS_Store" -type f -delete

deep-clean: clean

	@echo "    --------------------------------------------------"
	@echo "    Deep cleaning ..."
	@echo "    --------------------------------------------------"
	@$(MAKE) --no-print-directory -C $(BASH) deep-clean || true
	@$(MAKE) --no-print-directory -C $(PYTHON) deep-clean || true
	@$(MAKE) --no-print-directory -C $(R) deep-clean || true
	@$(MAKE) --no-print-directory -C $(SIMULATION) deep-clean || true
	@$(MAKE) --no-print-directory -C $(MARKDOWN) deep-clean || true

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    install          Build deal.II and ASPECT"
	@echo "    uninstall        Uninstall deal.II and ASPECT"
	@echo "    2d-shell-models  Run ASPECT shell models"
	@echo "    2d-box-models    Run ASPECT box models"
	@echo "    adiabat          Generate adiabatic profiles for ASPECT models"
	@echo "    visualize        Visualize all results"
	@echo "    environments     Create Conda environments"
	@echo "    clean            Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean       Deep clean figures and results (use with caution!)"
	@echo "    help             Show this help message"
