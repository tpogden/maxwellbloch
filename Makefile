.DEFAULT_GOAL := all

# Tests -----------------------------------------------------------------------

test:
	uv run pytest -n auto

test_cov:
	uv run pytest --cov -n auto

# Lint / Format ---------------------------------------------------------------

lint:
	uv run ruff check .

format:
	uv run ruff format .

format_check:
	uv run ruff format --check .

# Docs ------------------------------------------------------------------------

docs_html:
	uv run sphinx-build docs docs/_build -b html

# Dist ------------------------------------------------------------------------

dist:
	uv build

.PHONY: dist

# Deploy ----------------------------------------------------------------------

deploy_pypi_test: clean_dist dist
	uv run twine upload --repository-url https://test.pypi.org/legacy/ dist/*

deploy_pypi_prod: clean_dist dist
	uv run twine upload dist/*

# All -------------------------------------------------------------------------

all: test_cov docs_html dist

# Clean -----------------------------------------------------------------------

QU_FILES = $(shell find . -type f -name '**.qu')

clean_qu:
	@echo 'Deleting all **.qu files...'
	@echo $(QU_FILES)
	rm $(QU_FILES)

clean_docs:
	rm -rf docs/_build

clean_dist:
	rm -rf dist/*
