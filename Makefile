.DEFAULT_GOAL := all

# Tests -----------------------------------------------------------------------

test:
	uv run pytest -n auto

test_cov:
	uv run pytest --cov -n auto

bench:
	uv run pytest tests/bench_propagation.py --benchmark-only --benchmark-autosave --verbose

bench_max: 
	uv run pytest tests/bench_max.py --benchmark-only  --benchmark-autosave --verbose

# Lint / Format ---------------------------------------------------------------

lint:
	uv run ruff check .

format:
	uv run ruff format .

format_check:
	uv run ruff format --check .

# Docs ------------------------------------------------------------------------

# Quarto renders from stored notebook outputs; use docs_execute to re-run them.
docs: docs_html

docs_api:
	cd docs && uv run quartodoc build

docs_html: docs_api
	cd docs && uv run quarto render

# Re-execute every notebook (slow) and store the outputs in the .ipynb files.
docs_execute: docs_api
	cd docs && uv run ./execute.sh

docs_serve:
	cd docs && uv run quarto preview

serve: docs_serve

# Dist ------------------------------------------------------------------------

dist:
	uv build

.PHONY: dist docs docs_api docs_html docs_execute docs_serve serve clean_docs

# Release (bump version, commit, tag — then push to trigger CI publish) ------
# Usage: make bump_patch / bump_minor / bump_major

bump_patch:
	uv run bump-my-version bump patch

bump_minor:
	uv run bump-my-version bump minor

bump_major:
	uv run bump-my-version bump major

# All -------------------------------------------------------------------------

all: test_cov docs_html dist

# Clean -----------------------------------------------------------------------

QU_FILES = $(shell find . -type f -name '**.qu')

clean_qu:
	@echo 'Deleting all **.qu files...'
	@echo $(QU_FILES)
	rm $(QU_FILES)

clean_docs:
	rm -rf docs/_site

clean_dist:
	rm -rf dist/*

clean: clean_docs clean_qu clean_dist
