# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
R := $(PROJECT_ROOT)/R
PYTHON := $(PROJECT_ROOT)/python
FIGS := $(PROJECT_ROOT)/figs
MARKDOWN := $(PROJECT_ROOT)/markdown

# Manuscript
MANUSCRIPT := $(MARKDOWN)/manuscript
TEMPLATE := $(MANUSCRIPT)/eisvogel.tex
INCLUDE_FILES_LUA := $(MANUSCRIPT)/include-files.lua
CSL := $(MANUSCRIPT)/agu.csl
BIB := $(MANUSCRIPT)/main.bib
DOC := $(MANUSCRIPT)/manuscript
PREV_DOC := $(MANUSCRIPT)/bak/manuscript-07-oct-2025
DIFF_DOC := $(MANUSCRIPT)/diff

# Cleanup targets
CLEAN := rheological-parameters.md $(DOC).tex $(DOC).pdf $(DIFF_DOC).tex $(DIFF_DOC).pdf
DEEP_CLEAN :=

# Targets
.PHONY: all diff manuscript update-date draw-figures write-tables clean deep-clean help

all: manuscript

diff: $(DIFF_DOC).pdf

$(DIFF_DOC).pdf: $(PREV_DOC).tex $(DOC).tex
	@echo "    --------------------------------------------------"
	@echo "    Rendering manuscript diff pdf with citations"
	@echo "    --------------------------------------------------"
	@echo " -> Comparing $(shell basename $(PREV_DOC).tex) with $(shell basename $(DOC).tex)"
	@latexdiff \
		--flatten \
		--math-markup=1 \
		--type=CFONT \
		--append-textcmd=caption \
		--add-to-config PICTUREENV=longtable \
		--graphics-markup=none \
		$(PREV_DOC).tex $(DOC).tex > $(DIFF_DOC).tex
	@echo " -> Compiling $(shell basename $(DIFF_DOC).tex) to PDF"
	@pdflatex -interaction=nonstopmode -shell-escape -output-directory $(MANUSCRIPT) $(DIFF_DOC).tex
	@pdflatex -interaction=nonstopmode -shell-escape -output-directory $(MANUSCRIPT) $(DIFF_DOC).tex
	@rm -f $(DIFF_DOC).aux $(DIFF_DOC).log $(DIFF_DOC).toc

manuscript: $(DOC).tex $(DOC).pdf

$(DOC).pdf: update-date
	@echo "    --------------------------------------------------"
	@echo "    Rendering manuscript pdf"
	@echo "    --------------------------------------------------"
	@echo " -> Rendering $(shell basename $(DOC).md) with $(shell basename $(TEMPLATE)) template"
	@pandoc "$(DOC).md" -o "$(DOC).pdf" --from markdown --template "$(TEMPLATE)" --syntax-highlighting=kate \
		--bibliography=$(BIB) --csl=$(CSL) --filter=pandoc-crossref --citeproc --lua-filter=$(INCLUDE_FILES_LUA) \
		--pdf-engine=pdflatex --pdf-engine-opt="-shell-escape" --pdf-engine-opt="-interaction=nonstopmode"

$(DOC).tex: update-date
	@echo "    --------------------------------------------------"
	@echo "    Rendering manuscript tex"
	@echo "    --------------------------------------------------"
	@echo " -> Rendering $(shell basename $(DOC).md) with $(shell basename $(TEMPLATE)) template"
	@pandoc "$(DOC).md" -s -o "$(DOC).tex" --from markdown --template "$(TEMPLATE)" --syntax-highlighting=kate \
		--bibliography=$(BIB) --csl=$(CSL) --filter=pandoc-crossref --citeproc --lua-filter=$(INCLUDE_FILES_LUA)

update-date:
	@echo " -> Updating date in $(shell basename $(DOC).md)"
	@sed -i '' -E "s/^date: \".*\"/date: \"$(shell date '+%d %B %Y')\"/" $(DOC).md

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
