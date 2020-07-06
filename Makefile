test_cov:
	pytest --cov -n auto

QU_FILES = $(shell find . -type f -name '**.qu')

clean_qu:
	@echo 'Deleting all **.qu files...'
	@echo $(QU_FILES)
	rm $(QU_FILES)

docs_html:
	sphinx-build docs docs/_build -b html
