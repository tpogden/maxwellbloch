test_cov:
	pytest --cov -n 4

# clean_qu:
# 	@echo 'Deleting all *.qu files...'
# 	find . -type f -name '*.qu' 
	#-delete

docs_html:
	sphinx-build docs docs/_build -b html
