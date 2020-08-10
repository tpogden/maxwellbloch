.DEFAULT_GOAL := all

# Tests -----------------------------------------------------------------------

test:
	pytest -n auto

test_cov:
	pytest --cov -n auto

# Docs ------------------------------------------------------------------------

docs_html:
	sphinx-build docs docs/_build -b html

# Dist ------------------------------------------------------------------------

dist:
	python setup.py sdist --formats=gztar bdist_wheel

.PHONY: dist

# Deploy ----------------------------------------------------------------------

deploy_pypi_test: clean_dist dist
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

deploy_pypi_prod: clean_dist dist
	twine upload dist/*

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
