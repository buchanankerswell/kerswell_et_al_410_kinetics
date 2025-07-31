# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
R := $(PROJECT_ROOT)/R
PYTHON := $(PROJECT_ROOT)/python
FIGS := $(PROJECT_ROOT)/figs
MARKDOWN := $(PROJECT_ROOT)/markdown

# Manuscript
TEMPLATE := $(MARKDOWN)/manuscript/eisvogel.tex
INCLUDE_FILES_LUA := $(MARKDOWN)/manuscript/include-files.lua
CSL := $(MARKDOWN)/manuscript/agu.csl
BIB := $(MARKDOWN)/manuscript/main.bib
DOC := $(MARKDOWN)/manuscript/manuscript

# Cleanup targets
CLEAN := rheological-parameters.md
DEEP_CLEAN :=

# Targets
.PHONY: all manuscript draw-figures write-tables clean deep-clean help

all: manuscript

manuscript: $(DOC).pdf

$(DOC).pdf:
	@echo "    --------------------------------------------------"
	@echo "    Rendering manuscript"
	@echo "    --------------------------------------------------"
	@echo " -> Updating date in $(shell basename $(DOC).md)"
	@sed -i '' -E "s/^date: \".*\"/date: \"$(shell date '+%d %B %Y')\"/" $(DOC).md
	@echo " -> Rendering $(shell basename $(DOC).md) with $(shell basename $(TEMPLATE)) template"
	@pandoc "$(DOC).md" -o "$(DOC).pdf" --from markdown --template "$(TEMPLATE)" --listings \
		--bibliography=$(BIB) --csl=$(CSL) --filter=pandoc-crossref --citeproc \
		--lua-filter=$(INCLUDE_FILES_LUA) --pdf-engine=pdflatex --pdf-engine-opt="-shell-escape" \
		--pdf-engine-opt="-interaction=nonstopmode"

draw-figures:
	@$(MAKE) --no-print-directory -C $(R) visualize
	@$(MAKE) --no-print-directory -C $(PYTHON) visualize

write-tables:
	@$(MAKE) --no-print-directory -C $(R) manuscript

clean:
	@for item in $(CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				case "$$ABS_PATH" in \
					$(PROJECT_ROOT)*) \
						echo " -> Safely removing\n   $$ABS_PATH"; \
						rm -rf "$$ABS_PATH";; \
					*) \
						echo " -- Skipping (outside project root)\n   $$ABS_PATH";; \
				esac; \
			else \
				echo " -- Skipping (not found)\n   $$1"; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

deep-clean:
	@for item in $(DEEP_CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				case "$$ABS_PATH" in \
					$(PROJECT_ROOT)*) \
						echo " -> Safely removing\n   $$ABS_PATH"; \
						rm -rf "$$ABS_PATH";; \
					*) \
						echo " -- Skipping (outside project root)\n   $$ABS_PATH";; \
				esac; \
			else \
				echo " -- Skipping (not found)\n   $$1"; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    clean       Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean  Deep clean figures and results (use with caution!)"
	@echo "    help        Show this help message"
