.DEFAULT_GOAL := all

# Tests -----------------------------------------------------------------------

test:
	uv run pytest -n auto

test_cov:
	uv run pytest --cov -n auto

bench:
	uv run pytest tests/bench_mb_solve.py --benchmark-only

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
	rm -rf docs/_build

clean_dist:
	rm -rf dist/*
