.PHONY: all run dashboard test

PYTHON := .venv/bin/python

all: run

run:
	$(PYTHON) -B run_pipeline.py data/inputs/$(NAME).aa $(ARGS)

dashboard:
	$(PYTHON) -B dash_app_integrado.py

test:
	$(PYTHON) -B -m pytest
